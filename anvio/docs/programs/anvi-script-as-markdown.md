A helper script to format TAB-delimited files in markdown for bloggers, tutorial writers, those who wish to share example anvi'o outputs as text on GitHub issues, and so on.

Anvi'o programs often generate TAB-delimited files. While this simple format is useful to pass around to other software or share with others, it is not easily interpretable in visual media. The purpose of this script is to make the sharing part simpler for platforms that can render markdown.

You can pipe any TAB-delimited content to this script:

```
cat file.txt | anvi-sript-as-markdown
```

## Examples

Assume a TAB-delmited file with many lines:

``` bash
wc -l additional_view_data.txt

301 additional_view_data.txt
```

Contents of which look like this:

``` bash
head -n 10 additional_view_data.txt

contig	categorical_1	categorical_2	text_layer_01	numerical	bars_main!A	bars_main!B	bars_main!C
backrest	b	y	nmwje	2.78	278	23	1
backward	b	x	bqmyujr psrd doefhi	2.49	249	52	2
backwind	b	y	hkfer lchpmzix	2.69	269	32	3
backyard	b	x	advoe bfkyhmg	2.05	205	96	4
bacteria	b	x	lqmcwn hywco	2.63	263	38	5
bacterin	b		vxqdmn	2.98	298	3	6
baetylus	b	x	fkgpydi owgyhfx xwlpj	2.19	219	82	7
bagpiped	b	y	ijmnur	2.12	212	89	8
balconet	b	y	ecizgs	2.89	289	12	9
```

### Default run

{{ codestart }}
head -n 10 additional_view_data.txt | %(anvi-script-as-markdown)s
{{ codestop }}

which is rendered as,

|**contig**|**categorical_1**|**categorical_2**|**text_layer_01**|**numerical**|**bars_main!A**|**bars_main!B**|**bars_main!C**|
|:--|:--|:--|:--|:--|:--|:--|:--|
|backrest|b|y|nmwje|2.78|278|23|1|
|backward|b|x|bqmyujr psrd doefhi|2.49|249|52|2|
|backwind|b|y|hkfer lchpmzix|2.69|269|32|3|
|backyard|b|x|advoe bfkyhmg|2.05|205|96|4|
|bacteria|b|x|lqmcwn hywco|2.63|263|38|5|
|bacterin|b||vxqdmn|2.98|298|3|6|
|baetylus|b|x|fkgpydi owgyhfx xwlpj|2.19|219|82|7|
|bagpiped|b|y|ijmnur|2.12|212|89|8|
|balconet|b|y|ecizgs|2.89|289|12|9|

### Limit the number of lines shown

``` bash
cat additional_view_data.txt | anvi-script-as-markdown --max-num-lines-to-show 10
```

which is rendered as,

|**contig**|**categorical_1**|**categorical_2**|**text_layer_01**|**numerical**|**bars_main!A**|**bars_main!B**|**bars_main!C**|
|:--|:--|:--|:--|:--|:--|:--|:--|
|backrest|b|y|nmwje|2.78|278|23|1|
|backward|b|x|bqmyujr psrd doefhi|2.49|249|52|2|
|backwind|b|y|hkfer lchpmzix|2.69|269|32|3|
|backyard|b|x|advoe bfkyhmg|2.05|205|96|4|
|bacteria|b|x|lqmcwn hywco|2.63|263|38|5|
|bacterin|b||vxqdmn|2.98|298|3|6|
|baetylus|b|x|fkgpydi owgyhfx xwlpj|2.19|219|82|7|
|bagpiped|b|y|ijmnur|2.12|212|89|8|
|balconet|b|y|ecizgs|2.89|289|12|9|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

### Code columns

``` bash
cat additional_view_data.txt | anvi-script-as-markdown --max-num-lines-to-show 10 \
                                                       --code-column contig
```

which is rendered as,

|**contig**|**categorical_1**|**categorical_2**|**text_layer_01**|**numerical**|**bars_main!A**|**bars_main!B**|**bars_main!C**|
|:--|:--|:--|:--|:--|:--|:--|:--|
|`backrest`|b|y|nmwje|2.78|278|23|1|
|`backward`|b|x|bqmyujr psrd doefhi|2.49|249|52|2|
|`backwind`|b|y|hkfer lchpmzix|2.69|269|32|3|
|`backyard`|b|x|advoe bfkyhmg|2.05|205|96|4|
|`bacteria`|b|x|lqmcwn hywco|2.63|263|38|5|
|`bacterin`|b||vxqdmn|2.98|298|3|6|
|`baetylus`|b|x|fkgpydi owgyhfx xwlpj|2.19|219|82|7|
|`bagpiped`|b|y|ijmnur|2.12|212|89|8|
|`balconet`|b|y|ecizgs|2.89|289|12|9|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

### Exclude columns from the output

``` bash
cat additional_view_data.txt | anvi-script-as-markdown --max-num-lines-to-show 10 \
                                                       --code-column contig \
                                                       --exclude-columns 'bars_main!A,bars_main!B,bars_main!C'
```

which is rendered as,

|**contig**|**categorical_1**|**categorical_2**|**text_layer_01**|**numerical**|
|:--|:--|:--|:--|:--|
|`backrest`|b|y|nmwje|2.78|
|`backward`|b|x|bqmyujr psrd doefhi|2.49|
|`backwind`|b|y|hkfer lchpmzix|2.69|
|`backyard`|b|x|advoe bfkyhmg|2.05|
|`bacteria`|b|x|lqmcwn hywco|2.63|
|`bacterin`|b||vxqdmn|2.98|
|`baetylus`|b|x|fkgpydi owgyhfx xwlpj|2.19|
|`bagpiped`|b|y|ijmnur|2.12|
|`balconet`|b|y|ecizgs|2.89|
|(...)|(...)|(...)|(...)|(...)|
