This is a script that transposes tab-delimited files. That's it. 

It's helpful to get your inputs to line up with the types of inputs that anvi'o expects. Some programs have the `--transpose` flag, which will run this program for you, but some don't, and that's when you'll have to run it yourself. 

For example, anvi'o expects %(view-data)s to have each column representing a sample. If the file that you want to integrate into your anvi'o project has the samples as rows and the data attribute as the columns, then you'll need to %(anvi-script-transpose-matrix)s it. 

### An Example Run 

If you have an input ile `INPUT.txt` that looks like this: 

    1   2   3   
    4   5   6   
    7   8   9
    10  11  12
    
And you run this:

{{ codestart }}
anvi-script-transpose-matrix -o INPUT_transposed.txt \
                             -i INPUT.txt 
{{ codestop }}

You'll get a file called `INPUT_transposed.txt` that looks like 

    1   4   7   10
    2   5   8   11
    3   6   9   12
    

