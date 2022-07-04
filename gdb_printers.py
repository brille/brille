class BrillePrinter:
    def __init__(self, val):
        self.val = val

    def method(self, named, *args):
        arg_string = ",".join(f"{arg}" for arg in args) if len(args) else ""
        eval_string = f"(*({self.val.type}*)({self.val.address})).{named}({arg_string})"
        return gdb.parse_and_eval(eval_string)

    def wrap_to_string(self, string):
        return str(string)

    def _to_string(self):
        return ""

    def to_string(self):
        return self.wrap_to_string(self._to_string())


class BrilleLatticePrinter(BrillePrinter):
    def _to_string(self):
        return str(self.method('string_repr'))


class BrilleArrayPrinter(BrillePrinter):
    def _to_string(self):
        string = str(self.method('to_string')).replace(r"\n", "\n")
        return f"\n{string}"


class BrilleLatVecPrinter(BrillePrinter):
    def _to_string(self):
        array = BrilleArrayPrinter(self.val)
        lattice = BrilleLatticePrinter(self.val['lattice'])
        return f"\n{lattice._to_string()}\n{array._to_string()}"


def array_printer(val):
    if "Array" in str(val.type): return BrilleArrayPrinter(val)

def lattice_printer(val):
    if "Lattice" in  str(val.type): return BrilleLatticePrinter(val)

def latvec_printer(val):
    if "LatVeC" in str(val.type): return BrilleLatVecPrinter(val)

# gdb.pretty_printers.append(lattice_printer)
# gdb.pretty_printers.append(array_printer)
# gdb.pretty_printers.append(latvec_printer)