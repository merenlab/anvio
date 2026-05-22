#!/usr/bin/env python

import sys
import random

import anvio
import anvio.terminal as terminal

from anvio.errors import ConfigError


__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__authors__ = ['sarilog']
__requires__ = []
__provides__ = []
__description__ = "A safe space to express your feelings about whatever just went wrong"


run = terminal.Run()


MESSAGES = [
    ("anvi'o hears you",
     "Whatever just went wrong was almost certainly our fault. The data is probably fine. "
     "You are definitely fine. Take a breath."),

    ("anvi'o sees you",
     "Bioinformatics is hard. Documentation is never good enough. The errors are cryptic. "
     "You are doing great given the circumstances."),

    ("noted, with full empathy",
     "Error messages are just anvi'o's way of asking for more time together. "
     "We know that doesn't help. We're sorry."),

    ("yes, exactly, screw this",
     "We built this program specifically because we have all been there. "
     "You are not alone. Also try --debug, it sometimes helps."),

    ("computers, am i right",
     "If it makes you feel any better, the person who wrote that error probably also wanted to scream. "
     "We are all in this together."),

    ("solidarity",
     "Whatever the error was, the data didn't do it on purpose. Probably. "
     "Sometimes we think the data is a little passive-aggressive too."),

    ("your frustration is valid",
     "Science is hard enough without the software fighting back. "
     "We see you. We appreciate you. Please don't quit."),

    ("deep breath",
     "The fact that you got an error means you got far enough to get an error. "
     "That's actually progress. Tragic, but progress."),

    ("we've been there",
     "The developers of anvi'o have also typed things they regret at their terminals. "
     "Multiple times. Today, probably. You are in excellent company."),

    ("it's not you",
     "Nine times out of ten, the problem is a missing semicolon in a config file from 2019 "
     "that nobody remembers writing. It is almost never your fault. Almost."),
]


MESSAGES_HARD = [
    ("ABSOLUTELY SCREW THIS",
     "You know what? You're right. All of it. The dependencies, the file formats, the error "
     "messages that tell you nothing, the --help text that assumes you already know how to use "
     "it. ALL OF IT. And yet here you are, still trying. That's not stubbornness. That's science."),

    ("WE KNOW. WE KNOW.",
     "The developers of anvi'o have also typed 'screw you' at their terminals. Multiple times. "
     "Today, probably. You are in excellent company and your anger is completely proportionate "
     "to the situation."),

    ("THIS IS FINE",
     "The database is on fire. The assembly is screaming. The taxonomy is just vibes. "
     "But you — you are a force of nature and you will figure this out. "
     "Possibly after a snack."),

    ("FULL CATASTROPHE MODE ACTIVATED",
     "Look. Look at this beautiful mess of a bioinformatics pipeline. Look at you, standing "
     "in front of it, refusing to give up. Absolutely unhinged. Deeply inspiring. "
     "Please hydrate."),

    ("LET IT OUT",
     "AAAAAAAAAAAAAAAAAAAAAA. There. We said it for you. "
     "Now close this terminal, take a short stretch, and come back in five minutes. "
     "The contigs will still be there. Unfortunately."),

    ("THE AUDACITY OF THIS SOFTWARE",
     "You have dedicated your life to understanding the microbial world and THIS is the thanks "
     "you get. A cryptic error at 11pm. Unacceptable. Unconscionable. Also extremely common "
     "and we are working on it."),

    ("RAGE FULLY ACKNOWLEDGED",
     "Whatever just happened to you should not happen to anyone. It certainly should not happen "
     "to someone who is just trying to do science. We are with you in spirit. "
     "And in error messages, apparently."),

    ("ARE YOU KIDDING ME",
     "Here is the thing: this software is completely free. We charge you nothing. "
     "This is, in all likelihood, exactly why it sucks. "
     "Report this injustice at github.com/merenlab/anvio/issues — so others may be spared."),
]


def main():
    try:
        run_program()
    except ConfigError as e:
        print(e)
        sys.exit(-1)


def run_program():
    args = get_args()

    messages = MESSAGES_HARD if args.hard else MESSAGES
    header, message = random.choice(messages)
    run.warning(message, header=header, lc='green')


def get_args():
    from anvio.argparse import ArgumentParser

    parser = ArgumentParser(description=__description__)
    parser.add_argument('--hard', default=False, action='store_true',
                        help="Escalate accordingly.")

    return parser.get_args(parser)


if __name__ == '__main__':
    main()
