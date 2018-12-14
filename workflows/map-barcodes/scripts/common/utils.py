import itertools

# from https://stackoverflow.com/questions/3992735/python-generator-that-groups-another-iterable-into-groups-of-n
def chunks(iterable, n):
    iterable = iter(iterable)
    n_rest = n - 1
    for item in iterable:
        rest = itertools.islice(iterable, n_rest)
        yield itertools.chain((item,), rest)
