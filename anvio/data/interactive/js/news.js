/**
 * Javascript library to display anvi'o news.
 *
 *  Authors: Ozcan Esen
 *           Isaac Fink <iafink@uchicago.edu>
 *           A. Murat Eren <a.murat.eren@gmail.com>
 *
 * Copyright 2015-2021, The anvi'o project (http://anvio.org)
 *
 * Anvi'o is a free software. You can redistribute this program
 * and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * You should have received a copy of the GNU General Public License
 * along with anvi'o. If not, see <http://opensource.org/licenses/GPL-3.0>.
 *
 * @license GPL-3.0+ <http://opensource.org/licenses/GPL-3.0>
 */

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
