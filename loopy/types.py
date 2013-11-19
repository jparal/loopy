
class struct(dict):
    def __init__(self, *args, **kwargs):
        super(struct, self).__init__(*args, **kwargs)
        self.__dict__ = self
