$(document).ready(function() {
    if (getCookie('news_checked') == null) {
        checkNews();
    }
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
            createCookie('news_checked', 'yes', 1); // expiration is 1 days

            if (unread_count > 0) {
                $('#toggle-panel-right-3').css('color', '#FF0000');
            }
        }
    });
}

function newsMarkRead() {
    $('.blue-dot').remove();
    $('#toggle-panel-right-3').css('color', '#000000');
    createCookie('last_seen_hash', md5($('.news-item > h1')[0].textContent), -1);
}
