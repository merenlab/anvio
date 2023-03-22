# -*- coding: utf-8
# pylint: disable=line-too-long
"""Relations with the console output, Progress and Run classes"""

import os
import re
import sys
import time
import fcntl
import numpy as np
import struct
import pandas as pd
import termios
import datetime
import textwrap

from colored import fore, back, style
from collections import OrderedDict

import anvio
import anvio.dictio as dictio
import anvio.constants as constants

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


def pluralize(word, number, sfp="s", sfs=None, pfs=None, alt=None):
    """Pluralize a given word mindfully.

    We often run into a situation where the word of choice depends on the number of items
    it descripes and it can take a lot of extra space in the code. For instance, take this:

    >>> f"You have {num_sequences_in_fasta_file} sequences in your FASTA file."

    This will print "You have 1 sequences in your FASTA file" like an idiot when there is only
    one sequence. An alternative is to do something more elaborate:

    >>> f"You have {num_sequences_in_fasta_file} {'sequence' if num_sequences_in_fasta_file == 1 else 'sequences'}"

    Even though this will work beautifully, it works at the expense of the readability of the code
    for a minor inconvenience.

    THE PURPOSE of this function is to fix this problem in a more elegant fashion. The following call is
    equivalent to the second example:

    >>> f"You have {pluralize('sequence', num_sequences_in_fasta_file)} in your FASTA file."

    Alternatively, you can provide this function an `alt`, in which case it would return `word` for singular
    and `alt` for plural cases:

    >>> f"{pluralize('these do not', number, alt='this does not')}."

    Voila.


    Parameters
    ==========
    word: str
        The word to conditionally plurlize
    number: int
        The number of items the word intends to describe
    sfp: str, 's'
        Suffix for plural. The character that needs to be added to the end of the
        `word` if plural.
    sfs: str, None
        Suffix for singular. The same for `sfp` for singular case.
    pfs: str, None
        Prefix for singular. `pfs` will replace `1` in the final output (common
        parameters could be `pfs="a single" or pfs="only one"`).
    alternative: str, None
        If you provide an alternative, pluralize will discard every other parameter
        and will simply return `word` for singular, and `alt` for plural case.
    """

    plural = number != 1

    if plural:
        if alt:
            return alt
        else:
            return f"{pretty_print(number)} {word}{sfp}"
    else:
        if alt:
            return word
        else:
            if sfs:
                return f"{pretty_print(number)} {word}{sfs}"
            else:
                if pfs:
                    return f"{pfs} {word}"
                else:
                    return f"{number} {word}"


class Progress:
    def __init__(self, verbose=True):
        self.pid = None
        self.verbose = verbose
        self.terminal_width = None
        self.is_tty = sys.stdout.isatty()

        self.get_terminal_width()

        self.msg = None
        self.current = None

        self.progress_total_items = None
        self.progress_current_item = 0
        self.t = Timer(self.progress_total_items)

        self.LEN = lambda s: len(s.encode('utf-16-le')) // 2

        # if --no-progress or --quiet parameters were passed, OR, if we are not attached
        # to a real terminal, turn off progress outputs:
        if anvio.NO_PROGRESS or anvio.QUIET or not self.is_tty:
            self.verbose = False


    def get_terminal_width(self):
        try:
            self.terminal_width = max(get_terminal_size()[0], 60)
        except:
            # Getting the terminal size failed. It could be for many reasons: they may not have a
            # screen, they may be running TempleOS, etc. We respond by giving a generous terminal
            # width so that if they can see it at all, it truncates only the longest update messages.
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


    def update_pid(self, pid):
        self.pid = '%s %s' % (get_date(), pid)


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


    def update(self, msg, increment=False):
        self.msg = msg

        if not self.verbose:
            return

        if not self.pid:
            raise TerminalError('Progress with null pid will not update for msg "%s"' % msg)

        if increment:
            self.increment()

        self.clear()
        self.write('\r[%s] %s' % (self.pid, msg))


    def end(self, timing_filepath=None):
        """End the current progress

        Parameters
        ==========
        timing_filepath : str, None
            Store the timings of this progress to the filepath `timing_filepath`. File will only be
            made if a progress_total_items parameter was made during self.new()
        """

        if timing_filepath and self.progress_total_items is not None:
            self.t.gen_file_report(timing_filepath)

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

        self.single_line_prefixes = {0: '',
                                     1: '* ',
                                     2: '    - ',
                                     3: '        > '}

        if anvio.QUIET:
            self.verbose = False


    def log(self, line):
        if not self.log_file_path:
            self.warning("The run object got a logging request, but it was not inherited with "
                         "a log file path :(")
            return

        with open(self.log_file_path, "a") as log_file: log_file.write('[%s] %s\n' % (get_date(), CLEAR(line)))


    def write(self, line, quiet=False, overwrite_verbose=False):
        if self.log_file_path:
            self.log(line)

        if (self.verbose and not quiet) or (overwrite_verbose and not anvio.QUIET):
            try:
                sys.stderr.write(line)
            except:
                sys.stderr.write(line.encode('utf-8'))


    def info(self, key, value, quiet=False, display_only=False, overwrite_verbose=False, nl_before=0, nl_after=0, lc='cyan',
             mc='yellow', progress=None, align_long_values=True):
        """
        This function prints information nicely to the terminal in the form:
        key ........: value

        PARAMETERS
        ==========
        key : str
            what to print before the dots
        value : str
            what to print after the dots
        quiet : boolean
            the standard anvi'o quiet parameter which, if True, suppresses some output
        display_only : boolean
            if False, the key value pair is also stored in the info dictionary
        overwrite_verbose : boolean
            if True, downstream quiet parameters (though not the global --quiet) are ignored to produce more verbose output
        nl_before : int
            number of lines to print before the key-value line
        nl_after : int
            number of lines to print after the key-value line
        lc : color str
            the color of the label (key)
        mc : color str
            the color of the value
        progress : Progress instance
            provides the Progress bar to use
        align_long_values : boolean
            if True, values that are longer than the terminal width will be broken up into different lines that
            align nicely
        """
        if not display_only:
            self.info_dict[key] = value

        if value is None:
            value = "None"
        elif isinstance(value, bool) or isinstance(value, float) or isinstance(value, list):
            value = "%s" % value
        elif isinstance(value, str):
            value = remove_spaces(value)
        elif isinstance(value, int):
            value = pretty_print(value)

        label = constants.get_pretty_name(key)

        info_line = "%s%s %s: %s\n%s" % ('\n' * nl_before, c(label, lc),
                                         '.' * (self.width - len(label)),
                                         c(str(value), mc), '\n' * nl_after)
        if align_long_values:
            terminal_width = get_terminal_size()[0]
            wrap_width = terminal_width - self.width - 3
            wrapped_value_lines = textwrap.wrap(value, width=wrap_width, break_long_words=False, break_on_hyphens=False)
            if len(wrapped_value_lines) == 0:
                aligned_value_str = value
            else:
                aligned_value_str = wrapped_value_lines[0]
                for line in wrapped_value_lines[1:]:
                    aligned_value_str += "\n %s  %s" % (' ' * self.width, line)

            info_line = "%s%s %s: %s\n%s" % ('\n' * nl_before, c(label, lc),
                                             '.' * (self.width - len(label)),
                                             c(str(aligned_value_str), mc), '\n' * nl_after)

        if progress:
            progress.reset()
            self.write(info_line, overwrite_verbose=False, quiet=quiet)
            if progress.msg and progress.pid:
                progress.update(progress.msg)
        else:
            self.write(info_line, quiet=quiet, overwrite_verbose=overwrite_verbose)


    def info_single(self, message, overwrite_verbose=False, mc='yellow', nl_before=0, nl_after=0, cut_after=80, level=1, pretty_indentation=True, progress=None):
        if isinstance(message, str):
            message = remove_spaces(message)

        if level not in self.single_line_prefixes:
            raise TerminalError("the `info_single` function does not know how to deal with a level of %d :/" % level)

        if cut_after:
            if pretty_indentation:
                subsequent_indent = ''.join([' '] * len(self.single_line_prefixes[level]))
                message_line = c("%s%s\n" % (self.single_line_prefixes[level], textwrap.fill(str(message), cut_after, subsequent_indent=subsequent_indent)), mc)
            else:
                message_line = c("%s%s\n" % (self.single_line_prefixes[level], textwrap.fill(str(message), cut_after)), mc)
        else:
            message_line = c("%s%s\n" % (self.single_line_prefixes[level], str(message)), mc)

        message_line = ('\n' * nl_before) + message_line + ('\n' * nl_after)

        if progress:
            progress.reset()
            self.write(message_line, overwrite_verbose=False)
            if progress.msg and progress.pid:
                progress.update(progress.msg)
        else:
            self.write(message_line, overwrite_verbose=False)


    def warning(self, message, header='WARNING', lc='red', raw=False, overwrite_verbose=False, nl_before=0, nl_after=0, progress=None):
        if isinstance(message, str):
            message = remove_spaces(message)

        message_line = ''
        header_line = c("%s\n%s\n%s\n" % (('\n' * nl_before), header,
                                          '=' * (self.width + 2)), lc)
        if raw:
            message_line = c("%s\n\n%s" % ((message), '\n' * nl_after), lc)
        else:
            message_line = c("%s\n\n%s" % (textwrap.fill(str(message), 80), '\n' * nl_after), lc)

        if progress:
            progress.clear()
            self.write((header_line + message_line) if message else header_line, overwrite_verbose=overwrite_verbose)
            if progress.msg and progress.pid:
                progress.update(progress.msg)
        else:
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
    """Manages an ordered dictionary, where each key is a checkpoint name and value is a timestamp.

    Examples
    ========

    >>> from anvio.terminal import Timer
    >>> import time
    >>> t = Timer(); time.sleep(1)
    >>> t.make_checkpoint('checkpoint_name'); time.sleep(1)
    >>> timedelta = t.timedelta_to_checkpoint(timestamp=t.timestamp(), checkpoint_key='checkpoint_name')
    >>> print(t.format_time(timedelta, fmt = '{days} days, {hours} hours, {seconds} seconds', zero_padding=0))
    >>> print(t.time_elapsed())
    0 days, 0 hours, 1 seconds
    00:00:02

    >>> t = Timer(3) # 3 checkpoints expected until completion
    >>> for _ in range(3):
    >>>     time.sleep(1); t.make_checkpoint()
    >>>     print('complete: %s' % t.complete)
    >>>     print(t.eta(fmt='ETA: {seconds} seconds'))
    complete: False
    ETA: 02 seconds
    complete: False
    ETA: 01 seconds
    complete: True
    ETA: 00 seconds
    """
    def __init__(self, required_completion_score=None, initial_checkpoint_key=0, score=0):
        self.timer_start = self.timestamp()
        self.initial_checkpoint_key = initial_checkpoint_key
        self.last_checkpoint_key = self.initial_checkpoint_key
        self.checkpoints = OrderedDict([(initial_checkpoint_key, self.timer_start)])
        self.num_checkpoints = 0

        self.required_completion_score = required_completion_score
        self.score = score
        self.complete = False

        self.last_eta = None
        self.last_eta_timestamp = self.timer_start

        self.scores = {self.initial_checkpoint_key: self.score}


    def timestamp(self):
        return datetime.datetime.fromtimestamp(time.time())


    def timedelta_to_checkpoint(self, timestamp, checkpoint_key=None):
        if not checkpoint_key: checkpoint_key = self.initial_checkpoint_key
        timedelta = timestamp - self.checkpoints[checkpoint_key]
        return timedelta


    def make_checkpoint(self, checkpoint_key = None, increment_to = None):
        if not checkpoint_key:
            checkpoint_key = self.num_checkpoints + 1

        if checkpoint_key in self.checkpoints:
            raise TerminalError('Timer.make_checkpoint :: %s already exists as a checkpoint key. '
                                'All keys must be unique' % (str(checkpoint_key)))

        checkpoint = self.timestamp()

        self.checkpoints[checkpoint_key] = checkpoint
        self.last_checkpoint_key = checkpoint_key

        self.num_checkpoints += 1

        if increment_to:
            self.score = increment_to
        else:
            self.score += 1

        self.scores[checkpoint_key] = self.score

        if self.required_completion_score and self.score >= self.required_completion_score:
            self.complete = True

        return checkpoint


    def gen_report(self, title='Time Report', run=Run()):
        checkpoint_last = self.initial_checkpoint_key

        run.warning('', header=title, lc='yellow', nl_before=1, nl_after=0)

        for checkpoint_key, checkpoint in self.checkpoints.items():
            if checkpoint_key == self.initial_checkpoint_key:
                continue

            run.info(str(checkpoint_key), '+%s' % self.timedelta_to_checkpoint(checkpoint, checkpoint_key=checkpoint_last))
            checkpoint_last = checkpoint_key

        run.info('Total elapsed', '=%s' % self.timedelta_to_checkpoint(checkpoint, checkpoint_key=self.initial_checkpoint_key))


    def gen_dataframe_report(self):
        """Returns a dataframe"""

        d = {'key': [], 'time': [], 'score': []}
        for checkpoint_key, checkpoint in self.checkpoints.items():
            d['key'].append(checkpoint_key)
            d['time'].append(checkpoint)
            d['score'].append(self.scores[checkpoint_key])

        return pd.DataFrame(d)


    def gen_file_report(self, filepath):
        """Writes to filepath, will overwrite"""

        self.gen_dataframe_report().to_csv(filepath, sep='\t', index=False)


    def calculate_time_remaining(self, infinite_default = '∞:∞:∞'):
        if self.complete:
            return datetime.timedelta(seconds = 0)
        if not self.required_completion_score:
            return None
        if not self.score:
            return infinite_default

        time_elapsed = self.checkpoints[self.last_checkpoint_key] - self.checkpoints[0]
        fraction_completed = self.score / self.required_completion_score
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
        """Formats time

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
        """Run this and visually inspect its working"""

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
    """Time a block of code.

    This context manager times blocks of code, and calls run.info afterwards to report
    the time (unless quiet = True). See also time_program()

    Parameters
    ==========
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

    Examples
    ========

    >>> import time
    >>> import anvio.terminal as terminal
    >>> # EXAMPLE 1
    >>> with terminal.TimeCode() as t:
    >>>     time.sleep(5)
    ✓ Code finished successfully after 05s

    >>> # EXAMPLE 2
    >>> with terminal.TimeCode() as t:
    >>>     time.sleep(5)
    >>>     print(asdf) # undefined variable
    ✖ Code encountered error after 05s

    >>> # EXAMPLE 3
    >>> with terminal.TimeCode(quiet=True) as t:
    >>>     time.sleep(5)
    >>> print(t.time)
    0:00:05.000477
    """

    def __init__(self, success_msg=None, sc='green', fc='red', failure_msg=None, run=Run(), quiet=False, suppress_first=0):
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
    """A decorator used to time anvio programs.

    For a concrete example, see `bin/anvi-profile`.

    Examples
    ========

    >>> import anvio.terminal as terminal
    >>> @terminal.time_program
    >>> def main(args):
    >>>     <do stuff>
    >>> if __name__ == '__main__':
    >>>     <do stuff>
    >>>     main(args)
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


class TrackMemory(object):
    """Track the total memory over time

    Parameters
    ==========
    at_most_every : int or float, 5
        Memory is only calculated at most every 5 seconds, despite how many times self.measure is
        called
    """

    def __init__(self, at_most_every=5):
        self.t = None
        self.at_most_every = at_most_every


    def start(self):
        self.t = Timer(score=self._get_mem())
        return self.get_last(), self.get_last_diff()


    def measure(self):
        if self.t is None:
            raise TerminalError("TrackMemory :: You must start the tracker with self.start()")

        if self.t.timedelta_to_checkpoint(self.t.timestamp(), self.t.last_checkpoint_key) < datetime.timedelta(seconds = self.at_most_every):
            return False

        self.t.make_checkpoint(increment_to=self._get_mem())
        return True


    def gen_report(self):
        df = self.t.gen_dataframe_report().rename(columns={'score': 'bytes'}).set_index('key', drop=True)
        df['memory'] = df['bytes'].apply(self._format)
        return df


    def get_last(self):
        """Get the memory of the last measurement"""
        return self._format(self.t.scores[self.t.last_checkpoint_key])


    def get_last_diff(self):
        """Get the memory difference between the two latest measurements"""
        last_key = self.t.last_checkpoint_key

        if last_key == 0:
            return '+??'

        return self._format_diff(self._diff(last_key, last_key - 1))


    def _diff(self, key2, key1):
        return self.t.scores[key2] - self.t.scores[key1]


    def _format(self, mem):
        if np.isnan(mem):
            return '??'

        formatted = anvio.utils.human_readable_file_size(abs(mem))
        return ('-' if mem < 0 else '') + formatted


    def _format_diff(self, mem):
        if np.isnan(mem):
            return '+??'

        formatted = anvio.utils.human_readable_file_size(abs(mem))
        return ('-' if mem < 0 else '+') + formatted


    def _get_mem(self):
        mem = anvio.utils.get_total_memory_usage(keep_raw=True)

        if mem is None:
            return np.nan

        return mem


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


class Logger:
    """Utility class that makes it easier to use Anvio's nice logging in command runners."""
    def __init__(self, run=Run(), progress=Progress()):
        self.run = run
        self.progress = progress
