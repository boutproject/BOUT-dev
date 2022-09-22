class Enum(object):
    def __init__(self, name, vals, source):
        self.name = name.strip()
        self.vals = [x.strip() for x in vals.split(",")]
        if "None" in self.vals:
            self.vals.remove("None")
        self.lower = self.name.lower()
        self.extra = {}
        self.source = source.split("/")[-1]

    def __str__(self):
        return self.name


enums = {}
for inf in "../../../include/bout_types.hxx", "other_enums.hxx":
    with open(inf) as f:
        for line in f:
            line = line.strip()
            if line.startswith("enum"):
                assert line.startswith("enum class ")
                what = line[5 + 6 :]
                if __name__ == "__main__":
                    print(what, inf)
                name, vals = what.split("{")
                assert vals.endswith("};")
                enum = Enum(name, vals[:-2], inf)
                enums[enum.name] = enum
            if line.startswith("constexpr"):
                (enum, str), (val,) = [x.split() for x in line[10:].split("=")]
                val = val[:-1]  # remove ;
                val = val.split("::")[-1]
                if __name__ == "__main__":
                    print(enum, str, val)

enums["DIFF_METHOD"].extra = {
    "DIFF_DEFAULT": "deflt",
    "DEFAULT": "deflt",
    "DIFF_U1": "u1",
    "U1": "u1",
    "DIFF_U2": "u2",
    "U2": "u2",
    "DIFF_C2": "c2",
    "C2": "c2",
    "DIFF_W2": "w2",
    "W2": "w2",
    "DIFF_W3": "w3",
    "W3": "w3",
    "DIFF_C4": "c4",
    "C4": "c4",
    "DIFF_U3": "u3",
    "U3": "u3",
    "DIFF_FFT": "fft",
    "FFT": "fft",
    "DIFF_SPLIT": "split",
    "SPLIT": "split",
    "DIFF_S2": "s2",
    "S2": "s2",
}
enums["REGION"].extra = {
    "RGN_ALL": "all",
    "RGN_NOBNDRY": "nobndry",
    "RGN_NOX": "nox",
    "RGN_NOY": "noy",
    "RGN_NOZ": "noz",
}
enums["BRACKET_METHOD"].extra = {
    "BRACKET_STD": "standard",
    "STD": "standard",
    "BRACKET_SIMPLE": "simple",
    "SIMPLE": "simple",
    "BRACKET_ARAKAWA": "arakawa",
    "ARAKAWA": "arakawa",
    "BRACKET_CTU": "ctu",
    "CTU": "ctu",
    "BRACKET_ARAKAWA_OLD": "arakawa_old",
    "ARAKAWA_OLD": "arakawa_old",
}
enums["CELL_LOC"].extra = {
    "CELL_DEFAULT": "deflt",
    "DEFAULT": "deflt",
    "CELL_CENTRE": "centre",
    "CENTRE": "centre",
    "CELL_CENTER": "centre",
    "CENTER": "centre",
    "CELL_XLOW": "xlow",
    "XLOW": "xlow",
    "CELL_YLOW": "ylow",
    "YLOW": "ylow",
    "CELL_ZLOW": "zlow",
    "ZLOW": "zlow",
    "CELL_VSHIFT": "vshift",
    "VSHIFT": "vshift",
}

enums = list(enums.values())
