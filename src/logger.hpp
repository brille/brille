#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <filesystem>
#include <fstream>
#include <version>
#include <streambuf>
#include <memory>
#include <chrono>
#include <ctime>

#ifdef __cpp_lib_source_location
#include <source_location>
#define L_location location()
#else
#ifdef _MSC_VER
#define __PRETTY_FUNCTION__ __FUNCSIG__
#endif
#define L_location location(__FILE__,__LINE__, __PRETTY_FUNCTION__)
//#define L_location location(__FILE__,__LINE__, __func__)
#endif


namespace brille::logger{
enum Flags {
  clear,
  detab,
  timestamp,
  sync,
  endl,
  presync,
  silence,
  cout = 8,
  tabbed = 16,
  synced =32
};
enum Levels {
  hint, info, warn, error, critical
};

inline Flags operator += (Flags & l_flag, Flags r_flag){
  return l_flag = static_cast<Flags>(l_flag | r_flag);
}
inline Flags operator -= (Flags & l_flag, Flags r_flag){
  return l_flag = static_cast<Flags>(l_flag & ~r_flag);
}

#ifdef __cpp_lib_source_location
inline std::string location(const std::source_location& location = std::source_location::current()) {
  auto ss {std::stringstream {}};
  ss << location.file_name() << " : " << location.line() << " '" << location.function_name() << "'";
  return ss.str();
}
#else
inline std::string location(const char* file, int lineNo, const char* function) {
  auto ss = std::stringstream {};
  ss << file << " : " << lineNo << " '" << function << "'";
  return ss.str();
}
#endif

using Streamable = std::ostream;

class Logger {
public:
  void activate(bool makeActive=true){
    makeActive ? _flags -= silence : _flags += silence;
  }

  Flags addFlag(Flags flag) { return _flags += flag; }
  Flags removeFlag(Flags flag) { return _flags -= flag; }

  virtual void flush () {
    stream().flush();
    _flags -= presync;
  }

  virtual bool open() {return false;}

  template<class T> Logger& log(const T& value);

  Logger& operator << (Flags);
  Logger& operator << (decltype(std::endl<char, std::char_traits<char>>)) {
    return *this << endl;
  }
  Logger& operator << (decltype(std::hex) manip){
    stream() << manip;
    return *this;
  }
  Logger& operator << (decltype(std::setw) manip){
    stream() << manip;
    return *this;
  }

  virtual Streamable& stream();

  using ostreamPtr = Streamable*;

  virtual Logger* mirror_stream(ostreamPtr& mirrorStream) {
    mirrorStream = nullptr;
    return this;
  }

protected:
  explicit Logger(Flags initFlag = silence): _flags{initFlag} {}
  Logger(Flags initFlag, Streamable& = std::clog): _flags{initFlag} {}

  virtual Logger& logTime();

  template<class T> friend Logger& operator << (Logger& logger, T value);

  bool is_tabs() const {return _flags & tabbed || has_time();}
  bool is_null() const {return _flags == silence;}
  bool is_cout() const {return _flags & cout;}
  bool has_time() const {return (_flags & 7) == timestamp;}

  friend class FileNameGenerator;

  static tm* getTime();

  struct Log_date {
    unsigned char day;
    unsigned char month;
  } inline static log_date{0,0};

  Flags _flags = presync;
};

template<class T> Logger& Logger::log(const T& value){
  if (is_null()) return *this;
  auto streamPtr = &stream();
  Logger * logger = this;
  do {
    if (is_tabs()) *streamPtr << "\t";
    *streamPtr << value;
    logger = logger->mirror_stream(streamPtr);
  } while (streamPtr);
  removeFlag(timestamp);
  return *this;
}

class Null_Buff : public std::streambuf {
public:
  explicit Null_Buff() { setp(nullptr, nullptr); }
private:
  int_type overflow(int_type) override {
    return std::char_traits<char>::not_eof(0);
  }
} inline null_buff{};

inline Streamable null_ostream{&null_buff};

inline Streamable& Logger::stream() {return null_ostream;}

class ConsoleLogger: public Logger {
public:
  explicit ConsoleLogger(Flags init = silence, Streamable& ostream = std::clog):
    Logger{init, ostream}, _ostream{&ostream}{
    ostream.flush();
  }
  Streamable& stream() override {return is_null() ? Logger::stream() : *_ostream;}

  Logger* mirror_stream(ostreamPtr& mirrorStream) override {
    mirrorStream = mirrorStream == _ostream ? nullptr : _ostream;
    return this;
  }
protected:
  Streamable* _ostream{nullptr};
};

Logger& logger();

class FileNameGenerator{
public:
  static constexpr int FILE_NAME_LENGTH{8};
  FileNameGenerator(const std::filesystem::path& filepath);
  std::string stem() const {return _stem;}
  bool is_new_day(const Logger& l) const {return _day != l.log_date.day;}
  int day() const {return _day;}
  std::string operator()(const Logger& l);
private:
  std::string _stem;
  std::filesystem::path _path;
  unsigned char _day {0};
};


template<class MirrorBase = ConsoleLogger>
class FileLogger: public MirrorBase {
public:
  using path_t = std::filesystem::path;
  explicit FileLogger(const path_t& path): FileLogger{path, silence} {}
  FileLogger(const path_t& path, Flags init, Streamable& mirror = std::clog);
  FileLogger(const path_t& path, Flags init, Logger& mirror): FileLogger{path, init} {_mirror = &mirror;}
  Streamable& stream() override;
  void flush() override;
  Logger* mirror_stream(Logger::ostreamPtr& mirror) override;
  bool open() override;
private:
  Logger& logTime() override;
  FileNameGenerator _fileNameGenerator;
  Logger* _mirror{this};
  std::ofstream _file;
};
template<class MirrorBase> Streamable& FileLogger<MirrorBase>::stream() {
  if (MirrorBase::is_cout() || !open()){
    Logger::ostreamPtr streamPtr {&_file};
    mirror_stream(streamPtr);
    return *streamPtr;
  }
  return _file;
}
template<class MB> bool FileLogger<MB>::open(){
  if (_fileNameGenerator.is_new_day(*this)) _file.close();
  if (!_file.is_open()) _file.open(_fileNameGenerator(*this), std::ios::app); // append
  return _file.good();
}
template<class MB> Logger& FileLogger<MB>::logTime(){
  auto streamPtr {&stream()};
  auto* logger {mirror_stream(streamPtr)};
  while (streamPtr){
    *streamPtr << _fileNameGenerator.stem() << " ";
    logger = logger->mirror_stream(streamPtr);
  }
  MB::logTime();
  return *this;
}
template<class MB> void FileLogger<MB>::flush(){
  auto streamPtr {&stream()};
  auto* logger {mirror_stream(streamPtr)};
  while (streamPtr && logger != this){
    logger->flush();
    logger = logger->mirror_stream(streamPtr);
  }
  MB::flush();
  _file.flush();
}
template<class MB> Logger* FileLogger<MB>::mirror_stream(Logger::ostreamPtr& ms){
  bool isChained{this != _mirror};
  if (isChained){
    ms = &_mirror->stream();
    return _mirror;
  }
  return MB::mirror_stream(ms);
}

inline FileNameGenerator::FileNameGenerator(const std::filesystem::path &path): _path{path} {
  _stem = _path.filename().string();
  _stem.resize(FILE_NAME_LENGTH - 4);
  if (!path.has_extension()) _path += ".txt";
}
inline std::string FileNameGenerator::operator()(const Logger& l){
  if (l.log_date.day == 0) l.getTime();
  _day = l.log_date.day;
  auto fileName {std::stringstream{}};
  fileName << _stem << std::setfill('0') << std::setw(2)
           << (int)l.log_date.month << std::setw(2) << (int)_day;
  _path.replace_filename(fileName.str()) += _path.extension();
  return _path.string();
}

class RamBuffer: public std::streambuf{
public:
  RamBuffer(char* start, size_t size, Logger& l): _logger{&l} {setp(start, start+size);}
  void empty_buffer() {setp(pbase(), epptr());}
  auto start() const {return pbase();}
  auto pos() const {return pptr();}
private:
  int_type overflow(int_type ch) override {
    _logger->flush();
    sputc(static_cast<char>(ch));
    return std::char_traits<char>::not_eof(0);
  }
  Logger* _logger;
};

template<class MirrorBase=Logger>
class RamLogger: public FileLogger<MirrorBase> {
public:
  RamLogger(uint16_t size, const std::string& stem, Flags init, std::ostream& ostream = std::clog);
  std::ostream& stream() override {return _stream;}
  void flush() override;
private:
  std::unique_ptr<char[]> _ram_file;
  RamBuffer _ram_buffer;
  std::ostream _stream;
};
template<class MB> RamLogger<MB>::RamLogger(uint16_t size, const std::string& stem, Flags init, std::ostream& ostream)
: FileLogger<MB>{stem, init, ostream}, _ram_file{std::make_unique<char[]>(size)},
  _ram_buffer{_ram_file.get(), size, *this}, _stream{&_ram_buffer} {}

template<class MB> void RamLogger<MB>::flush(){
  for (char* c = _ram_buffer.start(); c < _ram_buffer.pos(); ++c)
    FileLogger<MB>::stream() << *c;
  _ram_buffer.empty_buffer();
}

}