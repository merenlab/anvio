$(document).ready(function() {
    checkNews();
});

function checkNews() {
    $('#news-panel-inner').empty();

    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/news',
        success: function(news) {
            var last_seen_hash = getCookie('last_seen_hash');
            var hash_found = false;
            var unread_count = 0;

            for (var i=0; i < news.length; i++) {
                var news_item = news[i];
                if (hash_found || last_seen_hash == md5(news_item['title'])) {
                    hash_found = true;
                } else {
                    unread_count++;
                }
                $('#news-panel-inner').append('<div class="news-item"> \
                                              <h1>' + ((hash_found) ? '' : '<span class="blue-dot">') + '</span>'+news_item['title']+'</h1> \
                                              <span class="news-date">'+news_item['date']+'</span>'+renderMarkdown(news_item['content'])+'</div>')
            }

            if (unread_count > 0) {
                $('#modNewsBadger').modal('show')
                $('#toggle-panel-right-3').css('color', '#FF0000');
                $('#toggle-panel-right-3').addClass('fading-button');

		$('#toggle-panel-right-3').mouseover(function() {
			if(!$('#news-panel').is(':visible')) {
				$(this).addClass('toggle-panel-right-pos-3');
                		toggleRightPanel('#news-panel');
                		$('#toggle-panel-right-3-inner').html('&#9658;');
                	}
		});
            } else {
		// remove upon starting interface in case 'last-seen-hash' cookie is not cleared
	    	$('#news-mark-read').remove();
            }
        }
    });
}

function newsMarkRead() {
    $('.blue-dot').remove();
    $('#toggle-panel-right-3').css('color', '#000000');
    $('#toggle-panel-right-3').removeClass('fading-button');
    $('#news-mark-read').remove();
    $('#toggle-panel-right-3').off('mouseover');
    createCookie('last_seen_hash', md5($('.news-item > h1')[0].textContent), -1);
}
