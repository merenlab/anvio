# -*- coding: utf-8
# pylint: disable=line-too-long
"""Relations with the console output, Progress and Run classes"""

import os
import re
import sys
import time
import fcntl
import struct
import termios
import datetime
import textwrap

from colored import fore, back, style
from collections import OrderedDict

import anvio.constants as constants
import anvio.dictio as dictio

from anvio.errors import TerminalError
from anvio.ttycolors import color_text as c

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "A. Murat Eren"
__email__ = "a.murat.eren@gmail.com"
__status__ = "Development"


# clean garbage garbage:
ansi_escape = re.compile(r'\x1B\[[0-?]*[ -/]*[@-~]')
non_ascii_escape = re.compile(r'[^\x00-\x7F]+')
CLEAR = lambda line: ansi_escape.sub('', non_ascii_escape.sub('', line.strip()))


class SuppressAllOutput(object):
    def __enter__(self):
        sys.stderr.flush()
        self.old_stderr = sys.stderr
        sys.stderr = open(os.devnull, 'w')
        sys.stdout.flush()
        self.old_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_value, traceback):
        sys.stderr.flush()
        sys.stderr = self.old_stderr
        sys.stdout.flush()
        sys.stdout = self.old_stdout


def remove_spaces(text):
    while True:
        if text.find("  ") > -1:
            text = text.replace("  ", " ")
        else:
            break

    return text


class Progress:
    def __init__(self, verbose=True):
        self.pid = None
        self.verbose = verbose
        self.terminal_width = None
        self.is_tty = sys.stdout.isatty()

        self.get_terminal_width()

        self.current = None

        self.progress_total_items = None
        self.progress_current_item = 0
        self.t = Timer(self.progress_total_items)

        self.LEN = lambda s: len(s.encode('utf-16-le')) // 2


    def get_terminal_width(self):
        # FIXME Program flow here is not clear. When does try fail?
        try:
            self.terminal_width = max(get_terminal_size()[0], 120)
        except:
            self.terminal_width = 120


    def new(self, pid, discard_previous_if_exists=False, progress_total_items=None):
        if self.pid:
            if discard_previous_if_exists:
                self.end()
            else:
                raise TerminalError("Progress.new() can't be called before ending the previous one (Existing: '%s', Competing: '%s')." % (self.pid, pid))

        if not self.verbose:
            return

        self.pid = '%s %s' % (get_date(), pid)
        self.get_terminal_width()
        self.current = None
        self.step = None
        self.progress_total_items = progress_total_items
        self.progress_current_item = 0
        self.t = Timer(self.progress_total_items)


    def increment(self, increment_to=None):
        if increment_to:
            self.progress_current_item = increment_to
        else:
            self.progress_current_item += 1

        self.t.make_checkpoint(increment_to = increment_to)


    def write(self, c, dont_update_current=False):
        eta_c = ' ETA: %s' % str(self.t.eta()) if self.progress_total_items else ''
        surpass = self.terminal_width - self.LEN(c) - self.LEN(eta_c)

        if surpass < 0:
            c = c[0:-(-surpass + 6)] + ' (...)'
        else:
            if not dont_update_current:
                self.current = c

            c += ' ' * surpass

        c += eta_c

        if self.verbose:
            if self.progress_total_items and self.is_tty:
                p_text = ''
                p_length = self.LEN(p_text)

                end_point = self.LEN(c) - self.LEN(eta_c)
                break_point = round(end_point * self.progress_current_item / self.progress_total_items)

                # see a full list of color codes: https://gitlab.com/dslackw/colored
                if p_length >= break_point:
                    sys.stderr.write(back.CYAN + fore.BLACK + c[:break_point] + \
                                     back.GREY_30 + fore.WHITE + c[break_point:end_point] + \
                                     back.CYAN + fore.CYAN + c[end_point] + \
                                     back.GREY_50 + fore.LIGHT_CYAN + c[end_point:] + \
                                     style.RESET)
                else:
                    sys.stderr.write(back.CYAN + fore.BLACK + c[:break_point - p_length] + \
                                     back.SALMON_1 + fore.BLACK + p_text + \
                                     back.GREY_30 + fore.WHITE + c[break_point:end_point] + \
                                     back.GREY_50 + fore.LIGHT_CYAN + c[end_point:] + \
                                     style.RESET)
                sys.stderr.flush()
            else:
                sys.stderr.write(back.CYAN + fore.BLACK + c + style.RESET)
                sys.stderr.flush()


    def reset(self):
        self.clear()


    def clear(self):
        if not self.verbose:
            return

        null = '\r' + ' ' * (self.terminal_width)
        sys.stderr.write(null)
        sys.stderr.write('\r')
        sys.stderr.flush()
        self.current = None
        self.step = None


    def append(self, msg):
        if not self.verbose:
            return
        self.write('%s%s' % (self.current, msg))


    def step_start(self, step, symbol="⚙ "):
        if not self.pid:
            raise TerminalError("You don't have an active progress to do it :/")

        if not self.current:
            raise TerminalError("You don't have a current progress bad :(")

        if self.step:
            raise TerminalError("You already have an unfinished step :( Here it is: '%s'." % self.step)

        if not self.verbose:
            return

        self.step = " / %s " % (step)

        self.write(self.current + self.step + symbol, dont_update_current=True)


    def step_end(self, symbol="👍"):
        if not self.step:
            raise TerminalError("You don't have an ongoing step :(")

        if not self.verbose:
            return

        self.write(self.current + self.step + symbol)

        self.step = None


    def update(self, msg):
        self.msg = msg

        if not self.verbose:
            return

        if not self.pid:
            raise TerminalError('Progress with null pid will not update for msg "%s"' % msg)

        self.clear()
        self.write('\r[%s] %s' % (self.pid, msg))


    def end(self):
        self.pid = None
        if not self.verbose:
            return
        self.clear()


class Run:
    def __init__(self, log_file_path=None, verbose=True, width=45):
        self.log_file_path = log_file_path

        self.info_dict = {}
        self.verbose = verbose
        self.width = width

        self.single_line_prefixes = {1: '* ',
                                     2: '    - ',
                                     3: '        > '}


    def log(self, line):
        if not self.log_file_path:
            self.warning("The run object got a logging request, but it was not inherited with\
                          a log file path :(")
            return

        with open(self.log_file_path, "a") as log_file: log_file.write('[%s] %s\n' % (get_date(), CLEAR(line)))


    def write(self, line, quiet=False, overwrite_verbose=False):
        if self.log_file_path:
            self.log(line)

        if (self.verbose and not quiet) or overwrite_verbose:
            try:
                sys.stderr.write(line)
            except:
                sys.stderr.write(line.encode('utf-8'))


    def info(self, key, value, quiet=False, display_only=False, nl_before=0, nl_after=0, lc='cyan', mc='yellow', progress=None):
        if not display_only:
            self.info_dict[key] = value

        if isinstance(value, bool):
            pass
        elif isinstance(value, str):
            value = remove_spaces(value)
        elif isinstance(value, int):
            value = pretty_print(value)

        label = constants.get_pretty_name(key)

        info_line = "%s%s %s: %s\n%s" % ('\n' * nl_before, c(label, lc),
                                         '.' * (self.width - len(label)),
                                         c(str(value), mc), '\n' * nl_after)

        if progress:
            progress.clear()
            self.write(info_line, quiet=quiet)
            progress.update(progress.msg)

        else:
            self.write(info_line, quiet=quiet)


    def info_single(self, message, mc='yellow', nl_before=0, nl_after=0, cut_after=80, level=1, progress=None):
        if isinstance(message, str):
            message = remove_spaces(message)

        if level not in self.single_line_prefixes:
            raise TerminalError("the `info_single` function does not know how to deal with a level of %d :/" % level)

        if cut_after:
            message_line = c("%s%s\n" % (self.single_line_prefixes[level], textwrap.fill(str(message), cut_after)), mc)
        else:
            message_line = c("%s%s\n" % (self.single_line_prefixes[level], str(message)), mc)

        message_line = ('\n' * nl_before) + message_line + ('\n' * nl_after)

        if progress:
            progress.clear()
            self.write(message_line)
            progress.update(progress.msg)

        else:
            self.write(message_line)


    def warning(self, message, header='WARNING', lc='red', raw=False, overwrite_verbose=False, nl_before=0, nl_after=0):
        if isinstance(message, str):
            message = remove_spaces(message)

        message_line = ''
        header_line = c("%s\n%s\n%s\n" % (('\n' * nl_before), header,
                                          '=' * (self.width + 2)), lc)
        if raw:
            message_line = c("%s\n\n%s" % ((message), '\n' * nl_after), lc)
        else:
            message_line = c("%s\n\n%s" % (textwrap.fill(str(message), 80), '\n' * nl_after), lc)

        self.write((header_line + message_line) if message else header_line, overwrite_verbose=overwrite_verbose)


    def store_info_dict(self, destination, strip_prefix=None):
        if strip_prefix:
            # mostly to get rid of output_dir prefix in output file names.
            # surprisingly enough, this is the best place to do it. live
            # and learn :/
            self.info_dict = dictio.strip_prefix_from_dict_values(self.info_dict, strip_prefix)

        dictio.write_serialized_object(self.info_dict, destination)


    def quit(self):
        if self.log_file_path:
            self.log('Bye.')


class Timer:
    """
    The premise of the class is to build an ordered dictionary, where each key is a checkpoint
    name and value is a timestamp.

    Examples
    ========

        from anvio.terminal import Timer
        import time
        t = Timer(); time.sleep(1)
        t.make_checkpoint('checkpoint_name'); time.sleep(1)
        timedelta = t.timedelta_to_checkpoint(timestamp=t.timestamp(), checkpoint_key='checkpoint_name')
        print(t.format_time(timedelta, fmt = '{days} days, {hours} hours, {seconds} seconds', zero_padding=0))
        print(t.time_elapsed())

        >>> 0 days, 0 hours, 1 seconds
        >>> 00:00:02

        t = Timer(3) # 3 checkpoints expected until completion
        for _ in range(3):
            time.sleep(1); t.make_checkpoint()
            print('complete: %s' % t.complete)
            print(t.eta(fmt='ETA: {seconds} seconds'))

        >>> complete: False
        >>> ETA: 02 seconds
        >>> complete: False
        >>> ETA: 01 seconds
        >>> complete: True
        >>> ETA: 00 seconds
    """
    def __init__(self, required_completion_score = None):
        self.timer_start = self.timestamp()
        self.last_checkpoint_key = 0
        self.checkpoints = OrderedDict([(self.last_checkpoint_key, self.timer_start)])
        self.num_checkpoints = 0

        self.required_completion_score = required_completion_score
        self.completion_score = 0
        self.complete = False

        self.last_eta = None
        self.last_eta_timestamp = self.timer_start


    def timestamp(self):
        return datetime.datetime.fromtimestamp(time.time())


    def timedelta_to_checkpoint(self, timestamp, checkpoint_key=0):
        timedelta = timestamp - self.checkpoints[checkpoint_key]
        return timedelta


    def make_checkpoint(self, checkpoint_key = None, increment_to = None):
        if not checkpoint_key:
            checkpoint_key = self.num_checkpoints + 1

        if checkpoint_key in self.checkpoints:
            raise TerminalError('Timer.make_checkpoint :: %s already exists as a checkpoint key.\
                                 All keys must be unique' % (str(checkpoint_key)))

        checkpoint = self.timestamp()

        self.checkpoints[checkpoint_key] = checkpoint
        self.last_checkpoint_key = checkpoint_key

        self.num_checkpoints += 1

        if increment_to:
            self.completion_score = increment_to
        else:
            self.completion_score += 1

        if self.required_completion_score and self.completion_score >= self.required_completion_score:
            self.complete = True

        return checkpoint


    def calculate_time_remaining(self, infinite_default = '∞:∞:∞'):
        if self.complete:
            return datetime.timedelta(seconds = 0)
        if not self.required_completion_score:
            return None
        if not self.completion_score:
            return infinite_default

        time_elapsed = self.checkpoints[self.last_checkpoint_key] - self.checkpoints[0]
        fraction_completed = self.completion_score / self.required_completion_score
        time_remaining_estimate = time_elapsed / fraction_completed - time_elapsed

        return time_remaining_estimate


    def eta(self, fmt=None, zero_padding=0):
        # Calling format_time hundreds or thousands of times per second is expensive. Therefore if
        # eta was called within the last half second, the previous ETA is returned without further
        # calculation.
        eta_timestamp = self.timestamp()
        if eta_timestamp - self.last_eta_timestamp < datetime.timedelta(seconds = 0.5) and self.num_checkpoints > 0:
            return self.last_eta

        eta = self.calculate_time_remaining()
        eta = self.format_time(eta, fmt, zero_padding) if isinstance(eta, datetime.timedelta) else str(eta)

        self.last_eta = eta
        self.last_eta_timestamp = eta_timestamp

        return eta


    def time_elapsed(self, fmt=None):
        return self.format_time(self.timedelta_to_checkpoint(self.timestamp(), checkpoint_key = 0), fmt=fmt)


    def format_time(self, timedelta, fmt = '{hours}:{minutes}:{seconds}', zero_padding = 2):
        """
            Examples of `fmt`. Suppose the timedelta is seconds = 1, minutes = 1, hours = 1.

            {hours}h {minutes}m {seconds}s  --> 01h 01m 01s
            {seconds} seconds               --> 3661 seconds
            {weeks} weeks {minutes} minutes --> 0 weeks 61 minutes
            {hours}h {seconds}s             --> 1h 61s
        """
        unit_hierarchy = ['seconds', 'minutes', 'hours', 'days', 'weeks']
        unit_denominations = {'weeks': 7, 'days': 24, 'hours': 60, 'minutes': 60, 'seconds': 1}

        if not fmt:
            # use the highest two non-zero units, e.g. if it is 7200s, use {hours}h{minutes}m
            seconds = int(timedelta.total_seconds())
            if seconds < 60:
                fmt = '{seconds}s'
            else:
                m = 1
                for i, unit in enumerate(unit_hierarchy):
                    if not seconds // (m * unit_denominations[unit]) >= 1:
                        fmt = '{%s}%s{%s}%s' % (unit_hierarchy[i-1],
                                                unit_hierarchy[i-1][0],
                                                unit_hierarchy[i-2],
                                                unit_hierarchy[i-2][0])
                        break
                    elif unit == unit_hierarchy[-1]:
                        fmt = '{%s}%s{%s}%s' % (unit_hierarchy[i],
                                                unit_hierarchy[i][0],
                                                unit_hierarchy[i-1],
                                                unit_hierarchy[i-1][0])
                        break
                    else:
                        m *= unit_denominations[unit]

        # parse units present in fmt
        format_order = []
        for i, x in enumerate(fmt):
            if x == '{':
                for j, k in enumerate(fmt[i:]):
                    if k == '}':
                        unit = fmt[i+1:i+j]
                        format_order.append(unit)
                        break

        if not format_order:
            raise TerminalError('Timer.format_time :: fmt = \'%s\' contains no time units.' % (fmt))

        for unit in format_order:
            if unit not in unit_hierarchy:
                raise TerminalError('Timer.format_time :: \'%s\' is not a valid unit. Use any of %s.'\
                                     % (unit, ', '.join(unit_hierarchy)))

        # calculate the value for each unit (e.g. 'seconds', 'days', etc) found in fmt
        format_values_dict = {}
        smallest_unit = unit_hierarchy[[unit in format_order for unit in unit_hierarchy].index(True)]
        units_less_than_or_equal_to_smallest_unit = unit_hierarchy[::-1][unit_hierarchy[::-1].index(smallest_unit):]
        seconds_in_base_unit = 1
        for a in [v for k, v in unit_denominations.items() if k in units_less_than_or_equal_to_smallest_unit]:
            seconds_in_base_unit *= a
        r = int(timedelta.total_seconds()) // seconds_in_base_unit

        for i, lower_unit in enumerate(unit_hierarchy):
            if lower_unit in format_order:
                m = 1
                for upper_unit in unit_hierarchy[i+1:]:
                    m *= unit_denominations[upper_unit]
                    if upper_unit in format_order:
                        format_values_dict[upper_unit], format_values_dict[lower_unit] = divmod(r, m)
                        break
                else:
                    format_values_dict[lower_unit] = r
                    break
                r = format_values_dict[upper_unit]

        format_values = [format_values_dict[unit] for unit in format_order]

        style_str = '0' + str(zero_padding) if zero_padding else ''
        for unit in format_order:
            fmt = fmt.replace('{%s}' % unit, '%' + '%s' % (style_str) + 'd')
        formatted_time = fmt % (*[format_value for format_value in format_values],)

        return formatted_time


    def _test_format_time(self):
        """
        Run this and visually inspect its working
        """
        run = Run()
        for exponent in range(1, 7):
            seconds = 10 ** exponent
            td = datetime.timedelta(seconds = seconds)

            run.warning('', header='TESTING %s' % td, lc='yellow')
            fmts = [
                None,
                "SECONDS {seconds}",
                "MINUTES {minutes}",
                "HOURS   {hours}",
                "DAYS    {days}",
                "WEEKS   {weeks}",

                "MINUTES {minutes} SECONDS {seconds}",
                "SECONDS {seconds} MINUTES {minutes}",
                "HOURS   {hours}   MINUTES {minutes}",
                "DAYS    {days}    HOURS   {hours}",
                "WEEKS   {weeks}   DAYS    {days}",
                "WEEKS   {weeks}   HOURS   {hours}",
                "WEEKS   {weeks}   MINUTES {minutes}",
                "DAYS    {days}    MINUTES {minutes}",
                "HOURS   {hours}   SECONDS {seconds}",

                "DAYS    {days}    MINUTES {minutes} SECONDS {seconds}",
                "WEEKS   {weeks}   HOURS {hours}     DAYS    {days}    SECONDS {seconds} MINUTES {minutes}",
            ]
            for fmt in fmts:
                run.info(str(fmt), self.format_time(td, fmt=fmt))


class TimeCode(object):
    """
    This context manager times blocks of code, and calls run.info afterwards to report
    the time (unless quiet = True). See also time_program()

    PARAMS
    ======
        sc: 'green'
            run info color with no runtime error
        success_msg: None
            If None, it is set to 'Code ran succesfully in'
        fc: 'green'
            run info color with runtime error
        failure_msg: None
            If None, it is set to 'Code failed within'
        run: Run()
            Provide a pre-existing Run instance if you want
        quiet: False,
            If True, run.info is not called and datetime object is stored
            as `time` (see examples)
        suppress_first: 0,
            Supress output if code finishes within this many seconds.

    EXAMPLES
    ========

        import time
        import anvio.terminal as terminal

        # EXAMPLE 1
        with terminal.TimeCode() as t:
            time.sleep(5)

        >>> ✓ Code finished successfully after 05s


        # EXAMPLE 2
        with terminal.TimeCode() as t:
            time.sleep(5)
            print(asdf) # undefined variable

        >>> ✖ Code encountered error after 05s

        # EXAMPLE 3
        with terminal.TimeCode(quiet=True) as t:
            time.sleep(5)
        print(t.time)

        >>> 0:00:05.000477
    """
    def __init__(self, sc='green', success_msg = None, fc='red', failure_msg = None, run = Run(), quiet = False, suppress_first = 0):
        self.run = run
        self.run.single_line_prefixes = {0: '✓ ', 1: '✖ '}

        self.quiet = quiet
        self.suppress_first = suppress_first
        self.sc, self.fc = sc, fc
        self.s_msg, self.f_msg = success_msg, failure_msg

        self.s_msg = self.s_msg if self.s_msg else 'Code finished after '
        self.f_msg = self.f_msg if self.f_msg else 'Code encountered error after '


    def __enter__(self):
        self.timer = Timer()
        return self


    def __exit__(self, exception_type, exception_value, traceback):
        self.time = self.timer.timedelta_to_checkpoint(self.timer.timestamp())

        if self.quiet or self.time <= datetime.timedelta(seconds=self.suppress_first):
            return

        return_code = 0 if exception_type is None else 1

        msg, color = (self.s_msg, self.sc) if not return_code else (self.f_msg, self.fc)
        self.run.info_single(msg + str(self.time), nl_before=1, mc=color, level=return_code)


def time_program(program_method):
    """
    A decorator used to time anvio programs. See below for example.
    For a concrete example, see `bin/anvi-profile`.

    EXAMPLE
    =======

    import anvio.terminal as terminal

    @terminal.time_program
    def main(args):
        <do stuff>

    if __name__ == '__main__':
        <do stuff>
        main(args)
    """
    import inspect
    program_name = os.path.basename(inspect.getfile(program_method))

    TimeCode_params = {
        'success_msg': '%s took ' % program_name,
        'failure_msg': '%s encountered an error after ' % program_name,
        'suppress_first': 3, # avoid clutter when program finishes or fails within 3 seconds
    }

    def wrapper(*args, **kwargs):
        with TimeCode(**TimeCode_params):
            program_method(*args, **kwargs)
    return wrapper


def pretty_print(n):
    """Pretty print function for very big integers"""
    if not isinstance(n, int):
        return n

    ret = []
    n = str(n)
    for i in range(len(n) - 1, -1, -1):
        ret.append(n[i])
        if (len(n) - i) % 3 == 0:
            ret.append(',')
    ret.reverse()
    return ''.join(ret[1:]) if ret[0] == ',' else ''.join(ret)


def tabulate(*args, **kwargs):
    """
    Uses the function `tabulate` in the `tabulate` module to tabulate data. This function behaves
    almost identically, but exists because currently multiline cells that have ANSI colors break the
    formatting of the table grid. These issues can be tracked to assess the status of this bug and
    whether or not it has been fixed:

    https://bitbucket.org/astanin/python-tabulate/issues/170/ansi-color-code-doesnt-work-with-linebreak
    https://bitbucket.org/astanin/python-tabulate/issues/176/ansi-color-codes-create-issues-with

    Until then, this overwrites a function in the module to preserve formatting when using multiline
    cells with ANSI color codes.
    """
    import tabulate

    def _align_column(strings, alignment, minwidth=0, has_invisible=True, enable_widechars=False, is_multiline=False):
        strings, padfn = tabulate._align_column_choose_padfn(strings, alignment, has_invisible)
        width_fn = tabulate._choose_width_fn(has_invisible, enable_widechars, is_multiline)
        s_widths = list(map(width_fn, strings))
        maxwidth = max(max(s_widths), minwidth)
        if is_multiline:
            if not enable_widechars and not has_invisible:
                padded_strings = [
                    "\n".join([padfn(maxwidth, s) for s in ms.splitlines()])
                    for ms in strings]
            else:
                lines = [line.splitlines() for line in strings]
                lines_pad = [[(s, maxwidth + len(s) - width_fn(s)) for s in group]
                             for group in lines]
                padded_strings = ["\n".join([padfn(w, s) for s, w in group])
                                  for group in lines_pad]
        else:
            if not enable_widechars and not has_invisible:
                padded_strings = [padfn(maxwidth, s) for s in strings]
            else:
                s_lens = list(map(len, strings))
                visible_widths = [maxwidth - (w - l) for w, l in zip(s_widths, s_lens)]
                padded_strings = [padfn(w, s) for s, w in zip(strings, visible_widths)]
        return padded_strings

    tabulate._align_column = _align_column
    return tabulate.tabulate(*args, **kwargs)


def get_date():
    return time.strftime("%d %b %y %H:%M:%S", time.localtime())


def get_terminal_size():
    """function was taken from http://stackoverflow.com/a/566752"""
    def ioctl_GWINSZ(fd):
        try:
            cr = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ,
        '1234'))
        except:
            return None
        return cr
    cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not cr:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            cr = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass
    if not cr:
        try:
            cr = (os.environ['LINES'], os.environ['COLUMNS'])
        except:
            cr = (25, 80)
    return int(cr[1]), int(cr[0])
