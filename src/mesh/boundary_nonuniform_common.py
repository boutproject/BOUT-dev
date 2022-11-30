orders = range(2, 5)
boundaries = ["Dirichlet", "Neumann", "Free"]


class maybeopen:
    def __init__(self, argv):
        if len(argv) > 1:
            self.fn = argv[1]
        else:
            self.fn = None

    def __enter__(self):
        if self.fn:
            self.fd = open(self.fn, "w")
        return self

    def __call__(self, *args, file=None):
        if file:
            print(*args, file=file)
        if self.fn:
            self.fd.write(" ".join(args) + "\n")
        else:
            print(*args)

    def __exit__(self, *args):
        if self.fn:
            self.fd.__exit__(*args)
