#include "logger.hpp"

using namespace brille::logger;

Logger & Logger::operator << (Flags flag) {
  if (is_null()) return *this;
  switch (flag) {
  case timestamp: logTime(); break;
  case sync:
    _flags = static_cast<Flags>(_flags & synced);
    *this << " |F|\n";
    flush();
    break;
  case endl: {
    if (_flags & synced) {*this << " |SYNC|";}
    else if (_flags == presync) {*this << " |PRESYNC|";}
    auto streamPtr = &stream();
    Logger * logger = this;
    do {
      *streamPtr << "\n";
      logger = logger->mirror_stream(streamPtr);
    } while (streamPtr);
    if (_flags & synced || _flags == presync) flush();
  }
    [[fallthrough]];
  case clear:
    if (_flags != presync) {_flags = static_cast<Flags>(_flags & synced);}
    break;
  case synced: _flags += synced; break;
  case detab: removeFlag(tabbed); break;
  default:
    addFlag(flag);
  }
  return *this;
}

tm* Logger::getTime() {
  std::time_t now = std::time(nullptr);
  auto localTime = std::localtime(&now);
  log_date.day = localTime->tm_mday;
  log_date.month = localTime->tm_mon + 1;
  return localTime;
}

Logger& Logger::logTime() {
  *this << std::put_time(getTime(), "%y-%m-%d %H:%M:%S");
  _flags += timestamp;
  return *this;
}

Logger& logger() {
  static ConsoleLogger std_log {};
  return std_log;
}