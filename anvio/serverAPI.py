# -*- coding: utf-8
"""A simple anvi'server API."""

import json
import requests

from anvio.terminal import Run
from anvio.errors import AnviServerError

run = Run()

class AnviServerAPI:
    def __init__(self, args):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        # server variables:
        self.user = A('user')
        self.password = A('password')

        # endpoint variables
        self.api_url = A('api_url')

        self.delete_if_exists = A('delete_if_exists')

        # project variables:
        self.tree = A('tree')
        self.items_order = A('items_order')
        self.view_data = A('view_data')
        self.additional_layers = A('additional_layers')
        self.fasta_file = A('fasta_file')
        self.layers_order_file = A('layers_order_file')
        self.layers_information_file = A('layers_information_file')

        self.project_name = A('project_name')
        self.description = A('description')
        self.state = A('state')
        self.bins = A('bins')
        self.bins_info = A('bins_info')

        self.cookie = None

    def login(self):
        r = self.request(path='/accounts/login/',
                         method='post',
                         data= {'username': self.user, 'password': self.password},
                         allow_redirects = False)

        if 'sessionid' in r.cookies:
            self.cookie = r.cookies['sessionid']
            run.info_single("Successfully logged into anvi-server.")
        else:
            raise Exception("Login failed, check your usename and password.")

    def get_url(self, path):
        return self.api_url + path

    def get_csrf_token(self):
        r = self.request(path='/accounts/login/', method='GET')

        if 'csrftoken' in r.cookies:
            return r.cookies['csrftoken']

        return ''

    def request(self, path = '', method = 'get', data = {}, files = {}, cookies = {}, allow_redirects = True):
        if self.cookie:
            cookies.update({'sessionid': self.cookie })

        if method == 'post':
            token = self.get_csrf_token()
            data.update({'csrfmiddlewaretoken': token })
            cookies.update({'csrftoken': token })

        try:
            return requests.request(method,
                                    self.get_url(path),
                                    data = data,
                                    files = files,
                                    cookies = cookies,
                                    allow_redirects = allow_redirects)
        except Exception as e:
            raise AnviServerError(str(e))

    def push(self):
        data = {'name': self.project_name}

        if self.description:
            data['description'] = open(self.description, 'r').read()
        else:
            data['description'] = ''

        if self.delete_if_exists:
            data['delete-if-exists'] = True

        files = {}

        if self.view_data:
            files['data.txt'] = open(self.view_data, 'r')
        if self.tree:
            files['tree.txt'] = open(self.tree, 'r')
        if self.items_order:
            files['items-order.txt'] = open(self.items_order, 'r')
        if self.fasta_file:
            files['fasta.fa'] = open(self.fasta_file, 'r')
        if self.layers_information_file:
            files['samples-info.txt'] = open(self.layers_information_file, 'r')
        if self.layers_order_file:
            files['samples-order.txt'] = open(self.layers_order_file, 'r')
        if self.additional_layers:
            files['additional-layers.txt'] = open(self.additional_layers, 'r')
        if self.state:
            files['state.json'] = open(self.state, 'r')
        if self.bins and self.bins_info:
            files['bins.txt'] = open(self.bins, 'r')
            files['bins-info.txt'] = open(self.bins_info, 'r')

        r = self.request(path='/projects/new',
                         method='post',
                         data=data,
                         files=files)

        response = json.loads(r.text)

        if response['status'] == 0:
            run.info_single('Project \'%s\' pushed successfully.' % self.project_name)
        else:
            raise AnviServerError(response['message'])
