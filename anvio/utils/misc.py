import time
import hashlib

from anvio.errors import ConfigError


def get_predicted_type_of_items_in_a_dict(d, key):
    """Gets a dictionary `d` and a `key` in it, and returns a type function.

    It is a bit counter intuitive. dictionary should look like this:

        d = {'x': {'key': item, (...)},
             'y': {'key': item, (...)},
             (...),
            }

    This is a shitty function, but there was a real need for it, so here we are :/
    """

    items = [x[key] for x in d.values()]

    if not items:
        # there is nothing to see here
        return None

    try:
        if(set(items) == set([None])):
            # all items is of type None.
            return None
    except TypeError:
        # this means we are working with an unhashable type.
        # it is either list or dict. we will go through items
        # and return the type of first item that is not None:
        for item in items:
            if item == None:
                continue
            else:
                return type(item)

        # the code should never come to this line since if everything
        # was None that would have been captured by the try block and the
        # exception would have never been thrown, but here is a final line
        # just to be sure we are not moving on with the rest of the code
        # if we entered into this block:
        return None

    # if we are here, it means not all items are None, and they are not of
    # unhashable types (so they must be atomic types such as int, float, or str)
    not_float = False
    for item in items:
        try:
            float(item or 0)
        except ValueError:
            not_float = True
            break

    if not_float:
        return str
    else:
        for item in items:
            try:
                if int(item or 0) == float(item or 0):
                    continue
                else:
                    return float
            except ValueError:
                return float

        return int



def HTMLColorToRGB(colorstring, scaled=True):
    """ convert #RRGGBB to an (R, G, B) tuple """
    colorstring = colorstring.strip()
    if colorstring[0] == '#': colorstring = colorstring[1:]
    if len(colorstring) != 6:
        raise ValueError("input #%s is not in #RRGGBB format" % colorstring)
    r, g, b = colorstring[:2], colorstring[2:4], colorstring[4:]
    r, g, b = [int(n, 16) for n in (r, g, b)]

    if scaled:
        return (r / 255.0, g / 255.0, b / 255.0)
    else:
        return (r, g, b)



def get_random_colors_dict(keys):
    # FIXME: someone's gotta implement this
    # keys   : set(1, 2, 3, ..)
    # returns: {1: '#ffffff', 2: '#888888', 3: '#222222', ...}
    return dict([(k, None) for k in keys])



def get_ordinal_from_integer(num):
    """append 'st', 'nd', or 'th' to integer to make categorical. num must be integer"""
    return'%d%s' % (num, {11:'th', 12:'th', 13:'th'}.get(num%100, {1:'st', 2:'nd', 3:'rd'}.get(num%10,'th')))



def get_f_string_evaluated_by_dict(f_string, d):
    """A function to evaluate the contents of an f-string given a dictionary.

    This simple function enables the following, even when the variables in the f-string
    are not defined in a given context, but appear as keys in a dictionary:

        >>> d = {'bar': 'apple', 'foo': 'pear', 'num': 5}
        >>> f_string = "{num}_{bar}_or_{foo}"
        >>> print(f"{get_f_string_evaluated_by_dict(f_string, d)}")
            "5_apple_or_pear"

    This functionality enables to receive a user-defined f-string from the commandline,
    and interpret it into a meaningful string using a dictionary. This is similar to the
    following use from earlier days of Python, but it doesn't bother the user to know
    about variable types and deal with an annoying syntax:

        >>> d = {'bar': 'apple', 'foo': 'pear', 'num': 5}
        >>> print("%(num)d_%(bar)s_or_%(foo)s" % d)
            "5_apple_or_pear"
    """

    stringlets = [p.split('}') for p in f_string.split('{')]

    if any([len(s) == 1 or len(s[0]) == 0 for s in stringlets[1:]]):
        raise ConfigError("Your f-string syntax is not working for anvi'o :/ Perhaps you "
                          "forgot to open or close a curly bracket?")

    unrecognized_vars = [s[0] for s in stringlets[1:] if s[0] not in d]
    if len(unrecognized_vars):
        raise ConfigError(f"Some of the variables in your f-string does not occur in the source "
                          f"dictionary :/ Here is the list of those that are not matching to anything: "
                          f"{', '.join(unrecognized_vars)}. In the meantime, these are the known keys: "
                          f"{', '.join(d.keys())}.")

    return stringlets[0][0] + ''.join([f"{d[s[0]]}{s[1]}" for s in stringlets[1:]])



def get_time_to_date(local_time, fmt='%Y-%m-%d %H:%M:%S'):
    try:
        local_time = float(local_time)
    except ValueError:
        raise ConfigError("utils::get_time_to_date is called with bad local_time.")

    return time.strftime(fmt, time.localtime(local_time))



def get_hash_for_list(l):
    return 'hash' + str(hashlib.sha224(''.join(sorted(list(l))).encode('utf-8')).hexdigest()[0:8])



def get_filtered_dict(input_dict, item, accepted_values_set):
    # removes any entry from d, where the value of the 'item' of items in d does not match
    # with 'accepted_values'
    if not isinstance(accepted_values_set, type(set([]))):
        raise ConfigError("get_filtered_dict: values must be type of set([]).")

    filtered_dict = {}

    for entry_id in input_dict:
        if input_dict[entry_id][item] not in accepted_values_set:
            continue
        else:
            filtered_dict[entry_id] = input_dict[entry_id]

    return filtered_dict

