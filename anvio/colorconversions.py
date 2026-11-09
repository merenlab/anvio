#!/usr/bin/env python
"""Conversions between color representations, and the mixing that depends on them.

The scope of this module is deliberately narrow: converting a color between representations, and the
two operations — blending and tinting — that are only correct once a color has been converted. It
holds no palettes, no choices about which color means what, and no knowledge of anvi'o data. Code
that decides which color to assign to something belongs where that decision is made, not here.

The reason the module exists at all is that a color hex code is not a measure of light. sRGB spends
most of its 256 steps per channel on dark tones, where the eye can tell them apart, and stretches
them out among bright ones, where it cannot: '#808080' emits about 21% of the light of '#ffffff'
rather than 50% of it. Averaging or interpolating the stored numbers therefore mixes encodings
rather than colors, and comes out too dark. Every function here that combines colors converts out of
that encoding first and back afterwards.
"""

import math
import matplotlib.colors as mcolors

from typing import Iterable, Tuple


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2026, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "Samuel Miller"
__email__ = "samuelmiller10@gmail.com"
__status__ = "Development"


def srgb_to_linear(channel: float) -> float:
    """Undo the sRGB transfer function, taking a stored channel value to light intensity."""
    return channel / 12.92 if channel <= 0.04045 else ((channel + 0.055) / 1.055) ** 2.4

def linear_to_srgb(channel: float) -> float:
    """Apply the sRGB transfer function, taking light intensity back to a stored channel value."""
    return 12.92 * channel if channel <= 0.0031308 else 1.055 * channel ** (1 / 2.4) - 0.055

def hex_to_oklab(hexcode: str) -> Tuple[float, float, float]:
    """
    Convert a color hex code to Oklab coordinates.

    Oklab (https://bottosson.github.io/posts/oklab/) is a perceptual color space in which the
    straight-line average of two colors looks like the mix a reader expects, and in which evenly
    spaced coordinates look evenly spaced. Neither is true of hex codes themselves, for the reason
    given in this module's docstring.

    Parameters
    ==========
    hexcode : str
        A color hex code, or any other Matplotlib color specification.

    Returns
    =======
    Tuple[float, float, float]
        The Oklab lightness and the two chromatic coordinates.
    """
    red, green, blue = (srgb_to_linear(channel) for channel in mcolors.to_rgb(hexcode))
    long_cone = 0.4122214708 * red + 0.5363325363 * green + 0.0514459929 * blue
    medium_cone = 0.2119034982 * red + 0.6806995451 * green + 0.1073969566 * blue
    short_cone = 0.0883024619 * red + 0.2817188376 * green + 0.6299787005 * blue
    # The cube roots are what make the space perceptually uniform, and are taken with 'copysign' so
    # that a slightly negative cone response, which rounding can produce for a color at the edge of
    # the sRGB gamut, does not raise on a fractional power.
    long_root = math.copysign(abs(long_cone) ** (1 / 3), long_cone)
    medium_root = math.copysign(abs(medium_cone) ** (1 / 3), medium_cone)
    short_root = math.copysign(abs(short_cone) ** (1 / 3), short_cone)
    return (
        0.2104542553 * long_root + 0.7936177850 * medium_root - 0.0040720468 * short_root,
        1.9779984951 * long_root - 2.4285922050 * medium_root + 0.4505937099 * short_root,
        0.0259040371 * long_root + 0.7827717662 * medium_root - 0.8086757660 * short_root
    )

def oklab_to_hex(oklab: Tuple[float, float, float]) -> str:
    """
    Convert Oklab coordinates back to a color hex code, clipping to the sRGB gamut.

    Parameters
    ==========
    oklab : Tuple[float, float, float]
        Oklab lightness and the two chromatic coordinates, as 'hex_to_oklab' returns.

    Returns
    =======
    str
        A lowercase six-digit hex code.
    """
    lightness, green_red, blue_yellow = oklab
    long_root = lightness + 0.3963377774 * green_red + 0.2158037573 * blue_yellow
    medium_root = lightness - 0.1055613458 * green_red - 0.0638541728 * blue_yellow
    short_root = lightness - 0.0894841775 * green_red - 1.2914855480 * blue_yellow
    long_cone = long_root ** 3
    medium_cone = medium_root ** 3
    short_cone = short_root ** 3
    red = 4.0767416621 * long_cone - 3.3077115913 * medium_cone + 0.2309699292 * short_cone
    green = -1.2684380046 * long_cone + 2.6097574011 * medium_cone - 0.3413193965 * short_cone
    blue = -0.0041960863 * long_cone - 0.7034186147 * medium_cone + 1.7076147010 * short_cone
    # The sRGB cube is not convex in Oklab, so a mix or a tint of colors that are themselves in
    # gamut can land outside it: tinting a saturated red asks for a linear red channel of up to
    # about 1.11. Clipping pulls each channel to the nearest value it can hold, which is what makes
    # this total — 'rgb2hex' raises on a channel outside 0-1 rather than clamping it itself. The
    # cost is that a clipped color is not exactly the one asked for, though the error stays far
    # below what a reader could see at the ramp lengths anvi'o draws.
    return mcolors.rgb2hex(
        tuple(min(1.0, max(0.0, linear_to_srgb(channel))) for channel in (red, green, blue))
    )

def blend_hexcodes(hexcodes: Iterable[str]) -> str:
    """
    Blend colors by averaging them perceptually.

    The result is the color between the given ones: blending green and blue gives the shade a reader
    would name as being both. Note that different sets of colors can of course blend to the same
    result, so a caller that needs its blends to be distinguishable has to check them.

    Parameters
    ==========
    hexcodes : Iterable[str]
        The colors to blend.

    Returns
    =======
    str
        A lowercase six-digit hex code.
    """
    coordinates = [hex_to_oklab(hexcode) for hexcode in hexcodes]
    return oklab_to_hex(tuple(sum(axis) / len(coordinates) for axis in zip(*coordinates)))


def tint_hexcode(hexcode: str, fraction: float) -> str:
    """
    Return a color a fraction of the way from white to the given color.

    This is how a ramp of one hue is built: evenly spaced fractions give evenly spaced tints, since
    the interpolation is perceptual. A fraction of 0 is white itself, which is worth knowing where
    white carries a meaning of its own.

    Parameters
    ==========
    hexcode : str
        The color the tint runs to.

    fraction : float
        Where between white (0.0) and that color (1.0) to land.

    Returns
    =======
    str
        A lowercase six-digit hex code.
    """
    white = hex_to_oklab('#ffffff')
    color = hex_to_oklab(hexcode)
    return oklab_to_hex(
        tuple(start + (end - start) * fraction for start, end in zip(white, color))
    )
