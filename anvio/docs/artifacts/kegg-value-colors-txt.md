This artifact is a TAB-delimited text file that gives a color of your choosing to each **categorical value** of a %(kegg-reaction-txt)s or %(kegg-compound-txt)s file drawn by %(anvi-draw-kegg-pathways)s.

A value column holds categorical values when the run gives `--reaction-value-type categorical` or `--compound-value-type categorical`. Such values are names, such as `up`, `down` and `same`. Normally anvi'o gives them its default colors in sorted order. This file says exactly which color each value gets instead, so that a color means the same thing across every figure in a paper. It also sets the order of the values.

Pass it with `--reaction-value-colors` for the reaction layer, and `--compound-value-colors` for the compound layer. The two layers are colored independently, so each can take its own file.

## File format

|Column|Required|Description|
|:--|:--|:--|
|*a value column*|YES|The **first** column, which can be headed anything you like, such as `value` or `status`. Each row names one categorical value.|
|`color`|YES|A six-digit color hex code, such as `#1f77b4`. Color *names* are not accepted, since the colors that pathway maps reserve for their own unhighlighted elements are compared as hex codes.|

Here is an example file:

|status|color|
|:--|:--|
|up|#d62728|
|down|#1f77b4|
|same|#bcbd22|

Every value of the layer needs a color. Anvi'o will not fill in the gaps with its default colors, since that would put colors you chose and colors chosen for you in one legend. Values in this file that the layer does *not* have are reported and ignored, so one file can serve several runs. Two values cannot share a color, since the maps and the legend could not tell them apart.

A value is a whole name, so a comma in it is just part of the name. This differs from a %(kegg-membership-colors-txt)s file, where a comma joins the names of a combination.

## The order of the rows

The order of the rows is the order of the values. It is used in three places.

- **The legend.** The legend lists the values in this order, with the first value at the top.
- **Drawing.** Where map elements overlap, the first value is drawn on top. `--reaction-reverse-overlay`/`--compound-reverse-overlay` puts the last value on top instead.
- **Choosing one value.** Where several values meet, a rule chooses one. The rules `plurality_order` and `order` use this order. `plurality_order` gives a tie to the value that comes first. `order` always takes the value that comes first, so the order can rank the values, for example to have `up` win over `same` in a map element whose KOs have both.

## When this file does not apply

This file needs categorical values. It is refused unless the layer's value type is `categorical`.

With a %(groups-txt)s, it is also refused unless the layer's sample summary gives a rule. That is `--reaction-sample-summary` for the reaction layer and `--compound-sample-summary` for the compound layer. Without a rule, each group's map counts the group's samples. No map then shows the values.
