#include <pybind11/pybind11.h>
#include "debug.hpp"

void wrap_debug(pybind11::module &m){
    using namespace pybind11::literals;
    using namespace brille;
    m.def("emit", [&](){return !printer.silenced();});
    m.def("emit", [&](const bool emt){
      bool orig = !printer.silenced();
      info_update_if( emt && orig, "Keeping printing on");
      info_update_if(!emt && orig, "Turning off printing");
      printer.silenced(!emt);
      info_update_if( emt && !orig, "Printing turned on");
      info_update_if(!emt && !orig, "Keeping printing off"); // never printed
      return !printer.silenced();
    }, "brille::printer emit status"_a);

    m.def("emit_datetime", [&](){return printer.datetime();});
    m.def("emit_datetime", [&](const bool emt){
      bool orig = printer.datetime();
      info_update_if( emt && orig, "Keeping datetime printing on");
      info_update_if(!emt && orig, "Turning off datetime printing");
      printer.datetime(emt);
      info_update_if( emt && !orig, "Datetime printing turned on");
      info_update_if(!emt && !orig, "Keeping datetime printing off");
      return printer.datetime();
    }, "brille::printer emit datetime status"_a);
}
