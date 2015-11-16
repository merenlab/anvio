# -*- coding: utf-8
"""
    configure sendmail
"""

import smtplib
from email.mime.text import MIMEText


__author__ = "Tobias Paczian"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = "1.0"
__maintainer__ = "Tobias Paczian"
__email__ = "tobiaspaczian@googlemail.com"
__status__ = "Development"

class mailsetup:
    def __init__(self):
        self.login = 'tobiaspaczian@googlemail.com'
        self.password = 'Alderaan42'


    def sendEmail(self, recipient, subject, messageText):

        adminEmail = 'tobiaspaczian@googlemail.com'
        adminPassword = 'Alderaan42'
        
        s = smtplib.SMTP_SSL('smtp.gmail.com', '465')
        #s.starttls()
        s.login(adminEmail, adminPassword)
        msg = MIMEText(messageText)
        msg['Subject'] = subject
        msg['From'] = adminEmail
        msg['To'] = recipient
        s.sendmail(adminEmail, [recipient], msg.as_string())
        s.quit()

        return True

