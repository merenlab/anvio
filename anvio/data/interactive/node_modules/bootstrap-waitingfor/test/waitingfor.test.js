var expect = chai.expect;

describe('Dialog testing', function () {

	var headerSelector = 'body > div.modal > div > div > div.modal-header > h3',
		modalBodySelector = 'body > div.modal > div > div > div.modal-body',
		modalSelector = 'body > div.modal';

	describe('Simple dialog', function () {

		it('displays a dialog with a progress bar', function (done) {
			waitingDialog.show();
			setTimeout(function () {
				expect($(headerSelector).text()).to.equal('Loading');
				waitingDialog.hide();
				done();
			}, 700);
		});

		it('displays a dialog with a custom message and a progress bar', function (done) {
			setTimeout(function () {
				waitingDialog.show('Custom message');
				setTimeout(function () {
					expect($(headerSelector).text()).to.equal('Custom message');
					waitingDialog.hide();
					done();
				}, 700);
			}, 300);
		});

		it('displays a dialog with custom settings and a progress bar', function (done) {
			setTimeout(function () {
				waitingDialog.show('Custom message 2', {
					dialogSize: 'sm', 
					progressType: 'warning'
				});
				setTimeout(function () {
					expect($(modalBodySelector).find('.progress-bar-warning').length).to.equal(1);
					expect($(headerSelector).text()).to.equal('Custom message 2');
					waitingDialog.hide();
					done();
				}, 700);
			}, 300);
		});

		it('displays a dialog with a special onHide event', function (done) {
			setTimeout(function () {
				waitingDialog.show('Custom message 3', {
					onHide: function () {
						expect($(modalSelector).is(':visible')).to.be.false;
						done();
					}
				});
				setTimeout(function () {
					expect($(headerSelector).text()).to.equal('Custom message 3');
					waitingDialog.hide();
				}, 700);
			}, 300);
		});

		it('displays a dialog with a special onShow event', function (done) {
			setTimeout(function () {
				waitingDialog.show('Custom message 4', {
					onShow: function () {
						expect($(modalSelector).is(':visible')).to.be.true;
					}
				});
				setTimeout(function () {
					expect($(headerSelector).text()).to.equal('Custom message 4');
					waitingDialog.hide();
					done();
				}, 700);
			}, 300);
		});

		it('testing "message" method', function (done) {
			setTimeout(function () {
				waitingDialog.show('Custom message 5', {
					onShow: function () {
						expect(waitingDialog.message()).to.be.equal('Custom message 5');
						waitingDialog.message('Changed message');
						expect(waitingDialog.message()).to.be.equal('Changed message');
					}
				});
				setTimeout(function () {
					waitingDialog.hide();
					done();
				}, 700);
			}, 300);
		});
		
		it('testing "message" method with custom headerSize', function (done) {
			setTimeout(function () {
				waitingDialog.show('Custom message 6', {
					headerSize: 1,
					onShow: function () {
						expect(waitingDialog.message()).to.be.equal('Custom message 6');
					}
				});
				setTimeout(function () {
					waitingDialog.hide();
					done();
				}, 700);
			}, 300);
		});

	});
});
