def deduplicate_tuple(t):
    """remove duplicate entries from the given tuple (starting from the left)"""
    seen_entries = []
    result = ()
    for entry in t:
        if entry in seen_entries:
            continue
        seen_entries.append(entry)
        result += (entry,)
    return result


def has_duplicates(t):
    """Return True, if given tuple t has duplicate entries, otherwise return False."""
    for d in _duplicates(t):
        return True
    return False


def signed_sort(t):
    """Sort the given tuple, returning the (-1)**n as well, where n is the number of transpositions needed for sorting"""
    if len(t) <= 1:
        return (1, t)

    min_pos = _position_of_minimal_element(t)
    (s, sl) = signed_sort(_remove_element(t, min_pos))
    return (s * ((-1) ** min_pos), (t[min_pos],) + sl)


def _remove_element(t, i):
    l = list(t)
    del l[i]
    return tuple(l)


def _position_of_minimal_element(t):
    if len(t) == 0:
        return None
    if len(t) == 1:
        return 0
    result = 0
    for i in range(len(t)):
        if t[i] < t[result]:
            result = i
    return result

def _duplicates(t):
    """generator yielding duplicate entries of given tuple (starting from the left)"""
    seen_entries = []
    for entry in t:
        if entry in seen_entries:
            yield entry
        seen_entries.append(entry)
