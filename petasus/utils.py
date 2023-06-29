"""Some small utility functions"""


def listify(obj):
    """Convert obj to a list, without splitting strings"""
    try:
        _ = iter(obj)
    except TypeError:
        obj = [obj]
    else:
        if isinstance(obj, str):
            obj = [obj]

    return list(obj)
