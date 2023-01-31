class polystarPrinter:
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


class polystarArrayPrinter(polystarPrinter):
    def _to_string(self):
        string = str(self.method('to_string')).replace(r"\n", "\n")
        return f"\n{string}"
def array_printer(val):
    if "Array" in str(val.type): return polystarArrayPrinter(val)

# gdb.pretty_printers.append(array_printer)
