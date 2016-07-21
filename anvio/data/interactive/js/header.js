function showHeader (currentPage, more) {
    var html = [];

    var menuitems = [
	{ "icon": "home", "title": "Server Home", "newpage": false, "href": "/" },
	{ "icon": "picture", "title": "Demo", "newpage": true, "href": "http://anvi-server.org/public/meren/interface_demo_VII" },
	{ "icon": "align-right", "title": "Interface and Data Types", "newpage": true, "href": "http://merenlab.org/2016/02/27/the-anvio-interactive-interface/" },
	{ "icon": "cog", "title": "Anvi\'o Project", "newpage": true, "href": "http://merenlab.org/projects/anvio/" },
    ];
    
    html.push('<nav class="navbar"><div class="container-fluid"><div class="navbar-header">');
    var menu = [];
    for (var i=0; i<menuitems.length; i++) {
	menu.push('<span class="glyphicon glyphicon-'+menuitems[i].icon+'"></span> <a href="'+menuitems[i].href+'"'+(menuitems[i].newpage ? ' target=_blank' : '')+'>'+menuitems[i].title+'</a>');
    }
    if (currentPage) {
	menu.push('<span class="glyphicon glyphicon-'+currentPage.icon+'"></span> '+currentPage.title);
    }
    html.push(menu.join(' / '));
    html.push('<div style="clear: both;font-size: 12px; padding-left: 30px; opacity: 0.8;" id="version"></div></div>')
    if (more) {
	html.push(more);
    }
    html.push('</div></nav>');
    
    document.getElementById('header').innerHTML = html.join("\n");
    $.getJSON('/version', function(data) {
	document.getElementById('version').innerHTML = "Running on anvi'o <b>v" + data.server + "</b> &amp; anvi'server DB <b>v"+ data.database + "</b>.";
    })
}
