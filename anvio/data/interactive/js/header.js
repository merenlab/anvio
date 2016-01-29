function showHeader (currentPage, more) {
    var html = [];

    var menuitems = [
	{ "icon": "picture", "title": "Gallery", "newpage": true, "href": "#" },
	{ "icon": "align-right", "title": "Data Types", "newpage": true, "href": "#" },
	{ "icon": "heart-empty", "title": "Anvi\'o Project", "newpage": true, "href": "http://merenlab.org/projects/anvio/" },
    ];

    menuitems.push(currentPage);
    
    html.push('<nav class="navbar"><div class="container-fluid"><div class="navbar-header">');
    var menu = [];
    for (var i=0; i<menuitems.length; i++) {
	menu.push('<span class="glyphicon glyphicon-'+menuitems[i].icon+'"></span> <a href="'+menuitems[i].href+'"'+(menuitems[i].newpage ? ' target=_blank' : '')+'>'+menuitems[i].title+'</a>');
    }
    if (currentPage) {
	menu.push('<span class="glyphicon glyphicon-'+currentPage.icon+'"></span> '+currentPage.title);
    }
    html.push(menu.join(' / '));
    if (more) {
	html.push(more);
    }
    html.push('</div></div></nav>');
    
    document.getElementById('header').innerHTML = html.join("\n");
}
