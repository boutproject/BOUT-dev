# what we need from bunch

class Bunch:
    def __init__(self, **dict):
        for k in dict:
            setattr(self, k, dict[k])
