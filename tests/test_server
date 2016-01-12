#!/usr/bin/env python
# -*- coding: utf-8
"""Test script of the server interface.

Performs unit tests of all available server functionality, excluding the analysis interface."""

__author__ = "Tobias Paczian"
__copyright__ = "Copyright 2015, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__maintainer__ = "Tobias Paczian"
__email__ = "tobiaspaczian@googlemail.com"
__status__ = "Development"

import unittest
from selenium import webdriver
from selenium.webdriver.common.alert import Alert
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

class AnvioTestCase(unittest.TestCase):

    EMAIL = "tobiaspaczian@googlemail.com"
    URL = "0.0.0.0:8080"

    def setUp(self):
        self.browser = webdriver.Chrome()
        self.addCleanup(self.browser.quit)

    def testRegisterValid(self):
        self.browser.get('http://'+self.URL)
        self.browser.find_element_by_link_text('Register').click()
        WebDriverWait(self.browser, 5).until(EC.title_is('anvio account registration'))
        self.browser.find_element_by_id('inputFirstname').send_keys('Test')
        self.browser.find_element_by_id('inputLastname').send_keys('User')
        self.browser.find_element_by_id('inputAffiliation').send_keys('anvio')
        self.browser.find_element_by_id('inputLogin').send_keys('testuser')
        self.browser.find_element_by_id('inputPassword').send_keys('test')
        self.browser.find_element_by_id('inputRepeatPassword').send_keys('test')
        self.browser.find_element_by_id('inputEmail').send_keys(self.EMAIL)
        self.browser.find_element_by_id('submit').click()
        WebDriverWait(self.browser, 5).until(EC.presence_of_element_located((By.XPATH, '//div[@class="alert alert-success col-sm-6"]')))
        self.assertTrue(self.browser.find_element_by_xpath('//div[@class="alert alert-success col-sm-6"]'))

    def testForgotPasswordInvalid(self):
        self.browser.get('http://'+self.URL)
        self.browser.find_element_by_link_text('forgot password?').click()
        WebDriverWait(self.browser, 5).until(EC.title_is('anvio forgot password'))
        self.browser.find_element_by_id('inputEmail').send_keys('invalid@email.com')
        self.browser.find_element_by_id('submit').click()
        WebDriverWait(self.browser, 5).until(EC.alert_is_present())
        self.assertEqual('Resetting password failed: No user has been found for email address "invalid@email.com"', Alert(self.browser).text)
        
    def testForgotPasswordValid(self):
        self.browser.get('http://'+self.URL)
        self.browser.find_element_by_link_text('forgot password?').click()
        WebDriverWait(self.browser, 5).until(EC.title_is('anvio forgot password'))
        self.browser.find_element_by_id('inputEmail').send_keys('tobi@email')
        self.browser.find_element_by_id('submit').click()
        WebDriverWait(self.browser, 5).until(EC.presence_of_element_located((By.XPATH, '//div[@class="alert alert-success col-sm-6"]')))
        self.assertTrue(self.browser.find_element_by_xpath('//div[@class="alert alert-success col-sm-6"]'))

    def testLogin(self):
        self.browser.get('http://'+self.URL)
        self.browser.find_element_by_id('login').send_keys('tobi')
        self.browser.find_element_by_id('password').send_keys('test')
        self.browser.find_element_by_id('loginButton').click()
        WebDriverWait(self.browser, 5).until(EC.title_is('anvio user home'))
        self.assertIn('anvio user home', self.browser.title)

    def testSelectProject(self):
        self.browser.get('http://'+self.URL)
        self.browser.find_element_by_id('login').send_keys('tobi')
        self.browser.find_element_by_id('password').send_keys('test')
        self.browser.find_element_by_id('loginButton').click() 
        WebDriverWait(self.browser, 5).until(EC.title_is('anvio user home'))
        self.browser.find_element_by_link_text('t_proj_01').click()
        WebDriverWait(self.browser, 5).until(EC.title_is('t_proj_01'))
        self.assertIn('t_proj_01', self.browser.title)

if __name__ == '__main__':
    unittest.main(verbosity=2)
