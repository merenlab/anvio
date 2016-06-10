# -*- coding: utf-8
# pylint: disable=line-too-long
"""A simple anvi'server API."""

import os
import json
import textwrap
import requests


class AnviServerError(Exception):
    def __init__(self, e=None):
        self.error_type = "Anvi'server Error:"
        self.e = e
        Exception.__init__(self)
        return

    def __str__(self):
        max_len = max([len(l) for l in textwrap.fill(self.e, 80).split('\n')])
        while True:
            if self.e.find("  ") > -1:
                self.e = self.e.replace("  ", " ")
            else:
                break

        error_lines = ['%s%s' % (l, ' ' * (max_len - len(l))) for l in textwrap.fill(self.e, 80).split('\n')]

        error_message = ['%s %s' % (self.error_type, error_lines[0])]
        for error_line in error_lines[1:]:
            error_message.append('%s%s' % (' ' * (len(self.error_type) + 2), error_line))

        return '\n\n' + '\n'.join(error_message) + '\n\n'

    def clear_text(self):
        return '%s: %s' % (self.error_type, self.e)



class AnviServerAPI:
    def __init__(self, args):
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

        # server variables:
        self.user = A('user')
        self.hostname = A('hostname')
        self.password = A('password')
        self.port_number = A('port_number')

        # project variables:
        self.tree = A('tree')
        self.view_data = A('view_data')
        self.additional_layer = A('additional_layers')
        self.fasta_file = A('fasta_file')
        self.samples_information = A('samples_information')
        self.samples_order = A('samples_order')

        self.project_name = A('project_name')
        self.state = A('state')

        self.version = 1
        self.protocol = 'http'

        self.token = None
        self.logged_in = False

        self.allowed_targets = set(['login', 'upload', 'project'])
        self.methods = {'POST': requests.post,
                        'DELETE': requests.delete,
                        'GET': requests.get}

        self.URL = lambda target: '%s://%s:%d/%s' % (self.protocol, self.hostname, self.port_number, self.is_target_proper(target))


    def is_target_proper(self, target):
        if target in self.allowed_targets:
            return target
        else:
            raise AnviServerError, "'%s' does not seem to be a target recognized by the API :/" % target


    def get_files_dict(self, files):
        if not files:
            return {}

        files_dict = {}

        for file_name in files:
            file_path = files[file_name]

            if not file_path:
                continue

            if not os.path.exists(file_path):
                raise AnviServerError, "File '%s' does not exist." % file_path

            files_dict[file_name] = open(file_path)

        return files_dict


    def request(self, target, data={}, files={}, method='POST', continue_on_error=False):
        if method not in self.methods:
            raise AnviServerError, "Unknown method: '%s'" % method

        url = self.URL(target)

        headers = {"X-Requested-With": "XMLHttpRequest", "User-Agent": "Anvi'o/v2"}

        cookies = dict(anvioSession=self.token) if self.token else dict()

        try:
            response_object = self.methods[method](url, data=data, headers=headers, cookies=cookies, files=files)
        except Exception as e:
            raise AnviServerError, "Something went wrong while trying to connect to the host %s. Here is a more\
                                detailed and uglier report: '''%s'''" % (self.hostname, e)

        if not response_object.status_code == requests.codes.ok:
            raise AnviServerError, "Somewhing went wrong with the server, and it returned an error code of %d.\
                                The reason for this error can be probably found in the server logs. If you\
                                would like to help, you can send an e-mail to us to let us know about htis\
                                error by telling us the details of what you were doing. Sorry!" \
                                                                            % response_object.status_code

        server_response = json.loads(response_object.text)

        if not continue_on_error and server_response['status'] == 'error':
            raise AnviServerError, server_response['message']

        return server_response


    def delete_project(self, project, continue_on_error=False):
        self.check_login()

        data = {'project': project}

        return self.request('project', data=data, method="DELETE", continue_on_error=continue_on_error)


    def login(self):
        data = {'login': self.user, 'password': self.password}

        response = self.request('login', data)

        if [t for t in ['token', 'firstname', 'lastname'] if t not in response['data']]:
            raise AnviServerError, 'The response from the server for a login request was not what the\
                                API was expecting to get. Something weird is happening.'

        self.token = response['data']['token']
        self.logged_in = True

        return response


    def check_login(self):
        if not self.logged_in:
            raise AnviServerError, "You are not logged in. Please make a login() call first."


    def push(self, delete_if_exists=False):
        self.check_login()

        if delete_if_exists:
            self.delete_project(self.project_name, continue_on_error=True)

        files = self.get_files_dict({'dataFile': self.view_data,
                                     'fastaFile': self.fasta_file,
                                     'treeFile': self.tree,
                                     'samplesOrderFile': self.samples_order,
                                     'samplesInformationFile': self.samples_information})

        data = {'title': self.project_name}

        return self.request('upload', data, files=files)


