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
import os

from selenium import webdriver
from selenium.webdriver.common.alert import Alert
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

class AnvioTestCase(unittest.TestCase):

    USER = 'testuser'
    EMAIL = "testuser@anvio.org"
    URL = "0.0.0.0:8080"
    BASEPATH = os.path.abspath("../sandbox/files_for_manual_interactive/") + '/'

    def setUp(self):
        self.browser = webdriver.Chrome()
        self.addCleanup(self.browser.quit)

    def login(self):
        self.browser.get('http://'+self.URL)
        self.browser.find_element_by_id('login').send_keys(self.USER)
        self.browser.find_element_by_id('password').send_keys('test')
        self.browser.find_element_by_id('loginButton').click()

    def loginAdmin(self):
        self.browser.get('http://'+self.URL)
        self.browser.find_element_by_id('login').send_keys('tobi')
        self.browser.find_element_by_id('password').send_keys('test')
        self.browser.find_element_by_id('loginButton').click()

    def test01RegisterNewUser(self):
        self.browser.get('http://'+self.URL)
        self.browser.find_element_by_link_text('Register').click()
        WebDriverWait(self.browser, 5).until(EC.title_is('anvio account registration'))
        self.browser.find_element_by_id('inputFirstname').send_keys('Test')
        self.browser.find_element_by_id('inputLastname').send_keys('User')
        self.browser.find_element_by_id('inputAffiliation').send_keys('anvio')
        self.browser.find_element_by_id('inputLogin').send_keys(self.USER+'2')
        self.browser.find_element_by_id('inputPassword').send_keys('test')
        self.browser.find_element_by_id('inputRepeatPassword').send_keys('test')
        self.browser.find_element_by_id('inputEmail').send_keys(self.EMAIL)
        self.browser.find_element_by_id('submit').click()
        WebDriverWait(self.browser, 5).until(EC.presence_of_element_located((By.XPATH, '//div[@class="alert alert-success col-sm-6"]')))
        self.assertTrue(self.browser.find_element_by_xpath('//div[@class="alert alert-success col-sm-6"]'))

    def test02Login(self):
        self.login()
        WebDriverWait(self.browser, 5).until(EC.title_is('anvio user home'))
        self.assertIn('anvio user home', self.browser.title)

    def test03ForgotPasswordInvalid(self):
        self.browser.get('http://'+self.URL)
        self.browser.find_element_by_link_text('forgot password?').click()
        WebDriverWait(self.browser, 5).until(EC.title_is('anvio forgot password'))
        self.browser.find_element_by_id('inputEmail').send_keys('invalid@email.com')
        self.browser.find_element_by_id('submit').click()
        WebDriverWait(self.browser, 5).until(EC.alert_is_present())
        self.assertEqual('Resetting password failed: No user has been found for email address "invalid@email.com"', Alert(self.browser).text)
        
    def test04ForgotPasswordValid(self):
        self.browser.get('http://'+self.URL)
        self.browser.find_element_by_link_text('forgot password?').click()
        WebDriverWait(self.browser, 5).until(EC.title_is('anvio forgot password'))
        self.browser.find_element_by_id('inputEmail').send_keys(self.USER+'@email')
        self.browser.find_element_by_id('submit').click()
        WebDriverWait(self.browser, 5).until(EC.presence_of_element_located((By.XPATH, '//div[@class="alert alert-success col-sm-6"]')))
        self.assertTrue(self.browser.find_element_by_xpath('//div[@class="alert alert-success col-sm-6"]'))

    def test05CreateNewProjectWithAllFiles(self):
        self.login()
        WebDriverWait(self.browser, 5).until(EC.visibility_of_element_located((By.XPATH, '//button[@title="upload data files"]')))

        self.browser.find_element_by_xpath('//button[@title="upload data files"]').click()
        WebDriverWait(self.browser, 5).until(EC.visibility_of_element_located((By.ID, 'uploadTitle')))
        self.browser.find_element_by_id('uploadTitle').send_keys('test_project_with_all_files')
        self.browser.find_element_by_id('uploadDescription').send_keys('description of test project')
        self.browser.find_element_by_id('treeFileSelect').send_keys(self.BASEPATH+'tree.txt')
        self.browser.find_element_by_id('fastaFileSelect').send_keys(self.BASEPATH+'fasta.fa')
        self.browser.find_element_by_id('dataFileSelect').send_keys(self.BASEPATH+'view_data.txt')
        self.browser.find_element_by_id('samplesOrderFileSelect').send_keys(self.BASEPATH+'samples-order.txt')
        self.browser.find_element_by_id('samplesInformationFileSelect').send_keys(self.BASEPATH+'samples-information.txt')
        self.browser.find_element_by_id('uploadFiles').click()
        WebDriverWait(self.browser, 5).until(EC.presence_of_element_located((By.LINK_TEXT, 'test_project_with_all_files')))
        self.assertTrue(self.browser.find_element_by_link_text('test_project_with_all_files'))

    def test06CreateNewProjectWithMinimalInput(self):
        self.login()
        WebDriverWait(self.browser, 5).until(EC.visibility_of_element_located((By.XPATH, '//button[@title="upload data files"]')))

        self.browser.find_element_by_xpath('//button[@title="upload data files"]').click()
        WebDriverWait(self.browser, 5).until(EC.visibility_of_element_located((By.ID, 'uploadTitle')))
        self.browser.find_element_by_id('uploadTitle').send_keys('test_project_minimal')
        self.browser.find_element_by_id('uploadDescription').send_keys('description of test project')
        self.browser.find_element_by_id('treeFileSelect').send_keys(self.BASEPATH+'tree.txt')
        self.browser.find_element_by_id('dataFileSelect').send_keys(self.BASEPATH+'view_data.txt')
        self.browser.find_element_by_id('uploadFiles').click()
        WebDriverWait(self.browser, 5).until(EC.presence_of_element_located((By.LINK_TEXT, 'test_project_minimal')))
        self.assertTrue(self.browser.find_element_by_link_text('test_project_minimal'))
 
    def test07SelectProject(self):
        self.login()
        WebDriverWait(self.browser, 5).until(EC.presence_of_element_located((By.LINK_TEXT, 'test_project_with_all_files')))
        self.browser.find_element_by_link_text('test_project_with_all_files').click()
        WebDriverWait(self.browser, 5).until(EC.title_is('test_project_with_all_files'))
        self.assertIn('test_project_with_all_files', self.browser.title)

    def test08ShareProject(self):
        self.login()
        WebDriverWait(self.browser, 5).until(EC.presence_of_element_located((By.XPATH, '//button[@title="share project"]')))
        self.browser.find_element_by_xpath('//button[@title="share project"]').click()
        WebDriverWait(self.browser, 5).until(EC.visibility_of_element_located((By.ID, 'projectName')))
        self.browser.find_element_by_id('projectName').send_keys('testShare')
        self.browser.find_element_by_id('shareProject').click()
        WebDriverWait(self.browser, 5).until(EC.visibility_of_element_located((By.ID, 'projectDescription')))
        self.assertTrue(self.browser.find_element_by_id('projectDescription'))

    def test09AddAdditionalDataFile(self):
        self.login()
        WebDriverWait(self.browser, 5).until(EC.presence_of_element_located((By.XPATH, '//button[@title="add data"]')))
        self.browser.find_element_by_xpath('//button[@title="add data"]').click()
        WebDriverWait(self.browser, 5).until(EC.visibility_of_element_located((By.ID, 'additionalFileType')))

        self.browser.find_element_by_id('additionalFileSelect').send_keys(self.BASEPATH+'additional_view_data.txt')
        self.browser.find_element_by_id('uploadAdditional').click()
        WebDriverWait(self.browser, 5).until(EC.visibility_of_element_located((By.XPATH, '//button[@title="project settings"]')))
        self.browser.find_element_by_xpath('//button[@title="project settings"]').click()
        WebDriverWait(self.browser, 5).until(EC.presence_of_element_located((By.ID, 'pfc3')))
        self.assertTrue(self.browser.find_element_by_id('pfc3'))

    def test10DeleteProject(self):
        self.login()
        WebDriverWait(self.browser, 5).until(EC.presence_of_element_located((By.ID, 'projectDelete0')))
        self.browser.find_element_by_id('projectDelete0').click()
        WebDriverWait(self.browser, 5).until(EC.alert_is_present())
        Alert(self.browser).accept()
        WebDriverWait(self.browser, 5).until(EC.visibility_of_element_located((By.XPATH, '//button[@title="upload data files"]')))
        self.assertFalse(len(self.browser.find_elements_by_link_text('test_proj_01')))

    def test11DeleteUser(self):
        self.loginAdmin()
        WebDriverWait(self.browser, 5).until(EC.presence_of_element_located((By.LINK_TEXT, 'admin page')))
        self.browser.find_element_by_link_text('admin page').click()
        WebDriverWait(self.browser, 5).until(EC.presence_of_element_located((By.XPATH, '//button[@onclick="showUserDetails(\'testuser2\');"]')))
        self.browser.find_element_by_xpath('//button[@onclick="showUserDetails(\'testuser2\');"]').click()
        WebDriverWait(self.browser, 5).until(EC.presence_of_element_located((By.XPATH, '//button[@onclick="deleteUser(\'testuser2\')"]')))
        self.browser.find_element_by_xpath('//button[@onclick="deleteUser(\'testuser2\')"]').click()
        WebDriverWait(self.browser, 5).until(EC.alert_is_present())
        Alert(self.browser).send_keys('CONFIRM')
        Alert(self.browser).accept()
        WebDriverWait(self.browser, 5).until(EC.presence_of_element_located((By.ID, 'usertable')))
        self.assertFalse(len(self.browser.find_elements_by_xpath('//button[@onclick="showUserDetails(\'testuser2\');"]')))

if __name__ == '__main__':
    unittest.main(verbosity=3)
