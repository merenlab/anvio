//--------------------------------------------------------------------------------------------------
// http://stackoverflow.com/questions/3019278/any-way-to-specify-the-base-of-math-log-in-javascript
function log10(val) {
  return Math.log(val) / Math.LN10;
}

//--------------------------------------------------------------------------------------------------
// http://stackoverflow.com/questions/1303646/check-whether-variable-is-number-or-string-in-javascript
function isNumber (o) {
  return ! isNaN (o-0);
}

//--------------------------------------------------------------------------------------------------
function ctype_alnum (str)
{
    return (str.match(/^[a-z0-9]+$/i) != null);
}

//--------------------------------------------------------------------------------------------------
function strip(html)
{
   var tmp = document.createElement("DIV");
   tmp.innerHTML = html;
   return tmp.textContent || tmp.innerText || "";
}

//--------------------------------------------------------------------------------------------------
function clearMinMax(selectbox) 
{
    var id = $(selectbox).attr('id').replace('normalization', '');

    $('#min' + id).val('0').prop('disabled', true);
    $('#max' + id).val('0').prop('disabled', true);     
}

//--------------------------------------------------------------------------------------------------
// source: https://gist.github.com/cjthompson/9140248
function readableNumber(num) {
    if(num == 0)
        return 0;
    var s = ['', 'K', 'M', 'G'];
    var e = Math.floor(Math.log(num) / Math.log(1000));
    return (num / Math.pow(1000, e)).toPrecision(3) + s[e];
}

//--------------------------------------------------------------------------------------------------
function linePath(p0, p1)
{
    var path = 'M ' + p0['x'] + ' ' + p0['y'] + ' ' + p1['x'] + ' ' + p1['y'];
    return path;
}

//--------------------------------------------------------------------------------------------------
function distance(p0, p1)
{
    return Math.sqrt(Math.pow(p1['x'] - p0['x'],2) + Math.pow(p1['y'] + p0['y'],2));
}

// http://stackoverflow.com/questions/498970/how-do-i-trim-a-string-in-javascript
if (!String.prototype.trim)
{
    String.prototype.trim=function(){return this.replace(/^\s+|\s+$/g, '');};
}

Math.toRadians = function(degrees) {
    return degrees * Math.PI / 180;
};

Math.toDegrees = function(radians) {
    return radians * 180 / Math.PI;
};


//--------------------------------------------------------------------------------------------------
//  Iterator
//--------------------------------------------------------------------------------------------------
function NodeIterator(root)
{
    this.root = root;
    this.cur = null;
    this.stack = [];
}

//--------------------------------------------------------------------------------------------------
NodeIterator.prototype.Begin = function() 
{
    this.cur = this.root;
    while (this.cur.child)
    {
        this.stack.push(this.cur);
        this.cur = this.cur.child;
    }
    return this.cur;
}

//--------------------------------------------------------------------------------------------------
NodeIterator.prototype.Next = function() 
{
    if (this.stack.length == 0)
    {
        this.cur = null;
    }
    else
    {
        if (this.cur.sibling)
        {
            var p = this.cur.sibling;
            while (p.child)
            {
                this.stack.push(p);
                p = p.child;
            }
            this.cur = p;
        }
        else
        {
            this.cur = this.stack.pop();
        }
    }
    return this.cur;
}

//--------------------------------------------------------------------------------------------------
PreorderIterator.prototype = new NodeIterator;

function PreorderIterator()
{
    NodeIterator.apply(this, arguments)
};

//--------------------------------------------------------------------------------------------------
PreorderIterator.prototype.Begin = function() 
{
    this.cur = this.root;
    return this.cur;
}

//--------------------------------------------------------------------------------------------------
PreorderIterator.prototype.Next = function() 
{
    if (this.cur.child)
    {
        this.stack.push(this.cur);
        this.cur = this.cur.child;
    }
    else
    {
        while (this.stack.length > 0 && this.cur.sibling == null)
        {
            this.cur = this.stack.pop();
        }
        if (this.stack.length == 0)
        {
            this.cur = null;
        }
        else
        {
            this.cur = this.cur.sibling;
        }
    }
    return this.cur;
}
