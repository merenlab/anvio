$(document).ready(function() {
    checkNews();
});

function checkNews() { // only render unread news items
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/news',
        success: function(news) {
            var last_seen_hash = getCookie('last_seen_hash');
            var hash_found = false;
            var unread_count = 0;

            // let test_news_item = {
            //     content: 'beep boop',
            //     date: '15.11.2021',
            //     title : 'this is my test news'
            // }
            // news.unshift(test_news_item)

            for (var i=0; i < news.length; i++) {
                var news_item = news[i];
                if (hash_found || last_seen_hash == md5(news_item['title'])) {
                    hash_found = true;
                } else {
                    unread_count++;
                    $('#modNewsBadger-inner').append('<div class="news-item"> \
                                                  <h1>' + ((hash_found) ? '' : '<span class="blue-dot">') + '</span>'+news_item['title']+'</h1> \
                                                  <span class="news-date">'+news_item['date']+'</span>'+renderMarkdown(news_item['content'])+'</div>')
                }
            }
            if (unread_count > 0) {
                $('#modNewsBadger').modal('show')
            }
        }
    });
}

function showAllNews(){ // render all news items regardless of cookie content - called from interactive hamburger dropdown.
    $('#modNewsBadger-inner').empty()
    console.log('got here')
    $.ajax({
        type: 'GET',
        cache: false,
        url: '/data/news',
        success: function(news) {
            for (var i=0; i < news.length; i++) {
                var news_item = news[i];
                $('#modNewsBadger-inner').append('<div class="news-item"> \
                                              <h1>' + '</span>'+news_item['title']+'</h1> \
                                              <span class="news-date">'+news_item['date']+'</span>'+renderMarkdown(news_item['content'])+'</div>')
            }
            $('#modNewsBadger').modal('show')
        }
    });
}

function newsMarkRead() {
    $('.blue-dot').remove();
    createCookie('last_seen_hash', md5($('.news-item > h1')[0].textContent), -1);
}
