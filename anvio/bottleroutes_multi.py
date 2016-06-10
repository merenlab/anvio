# -*- coding: utf-8
# pylint: disable=line-too-long
"""
    Multi-User routes for bottle web server.

    Functions defined here are used from wihtin programs such as
    anvi-interactive, or anvi-refine.
"""

import json
import os

from bottle import redirect, static_file

import anvio
import anvio.utils as utils
import anvio.dbops as dbops
import anvio.terminal as terminal

run = terminal.Run()

def set_default_headers(response):
    response.set_header('Content-Type', 'application/json')
    response.set_header('Pragma', 'no-cache')
    response.set_header('Cache-Control', 'no-cache, no-store, max-age=0, must-revalidate')
    response.set_header('Expires', 'Thu, 01 Dec 1994 16:00:00 GMT')


def server_version(request, userdb, response):
    set_default_headers(response)
    return json.dumps({"database": anvio.__users_db_version__, "server": anvio.__version__})


def get_user(request, userdb, response):
    # check if we have a cookie
    if request.get_cookie('anvioSession'):

        # we have a cookie, check if it is valid
        user = userdb.get_user_for_token(request.get_cookie('anvioSession'))
        if user['status'] == 'ok':
            return user
        else:
            return False

    else:
        return False

def get_user_by_token(request, userdb, response):
    set_default_headers(response)
    return json.dumps(userdb.get_user_for_token(request.forms.get('token'), request.forms.get('verbose')))


def impersonate(request, userdb, response):
    set_default_headers(response)
    retval = get_user(request, userdb, response)
    if retval:
        if retval['data']['clearance'] == 'admin':
            imperson = userdb.get_user_for_login(request.forms.get('login'))
            if imperson['status'] == 'ok':
                response.set_header('Set-Cookie', 'anvioSession=' + imperson['data']["token"] + '; path=/; max-age=' + str(60 * 60 * 24 * 14))
            return json.dumps(imperson)


def change_clearance(request, userdb, response):
    set_default_headers(response)
    return json.dumps(userdb.change_clearance(request.forms.get('user'), request.forms.get('clearance'), get_user(request, userdb, response)))


def request_account(request, userdb, response):
    set_default_headers(response)
    return json.dumps(userdb.create_user(request.forms.get('firstname'), request.forms.get('lastname'), request.forms.get('email'), request.forms.get('login'), request.forms.get('password'), request.forms.get('affiliation'), request.environ.get('REMOTE_ADDR'), ))


def accept_user(request, userdb, response):
    retval = userdb.accept_user(request.query.login, request.query.code)
    if retval['status'] == 'ok':
        redirect('/app/accountOK.html')
    else:
        redirect('/app/accountBAD.html')


def delete_user(request, userdb, response):
    set_default_headers(response)
    user = get_user(request, userdb, response)
    if user:
        return json.dumps(userdb.delete_user(request.forms.get('user'), user))
    else:
        return '{ "status": "error", "data": null, "message": "you must be logged in as an admin to delete a user" }'

def reset_password(request, userdb, response):
    set_default_headers(response)
    return json.dumps(userdb.get_user_for_email(request.forms.get('email')))


def check_availability(request, userdb, response):
    set_default_headers(response)
    if request.forms.get('email'):
        user = userdb.get_user_for_email(request.forms.get('email'))
        if user['status'] == 'ok':
            return '{ "status": "error", "data": "email", "message": "this email is already taken" }'
        else:
            return '{ "status": "ok", "message": null, "data": "email" }'
    elif request.forms.get('login'):
        user = userdb.get_user_for_login(request.forms.get('login'))
        if user['status'] == 'ok':
            return '{ "status": "error", "message": "this login is already taken", "data": "login" }'
        else:
            return '{ "status": "ok", "message": null, "data": "login" }'


def change_password(request, userdb, response):
    set_default_headers(response)
    user = userdb.get_user_for_login(request.forms.get('login'))
    if user['status'] == 'ok':
        return json.dumps(userdb.change_password(user, request.forms.get('password')))
    else:
        return json.dumps(user)


def login_to_app(request, userdb, response):
    set_default_headers(response)
    retval = userdb.login_user(request.forms.get('login'), request.forms.get('password'))
    if retval['status'] == 'ok':
        response.set_header('Set-Cookie', 'anvioSession=' + retval['data']["token"] + '; path=/; max-age=' + str(60 * 60 * 24 * 14))

    return json.dumps(retval)


def logout_from_app(request, userdb, response):
    return json.dumps(userdb.logout_user(request.forms.get('login')))


def set_view_cookie(request, userdb, response, login, name, private):
    # get the absolute path for static directory under anvio

    static_dir = os.path.join(os.path.dirname(utils.__file__), 'data/interactive')
    response = static_file('index.html', root=static_dir)

    if not userdb.view_exists(login, name):
        redirect('/app/error.html?err=404')
        return

    if private:
        if not request.query.code:
            redirect('/app/error.html?err=401')
            return
        else:
            response.set_header('Set-Cookie', 'anvioView=' + login + '|' + name + '|' + request.query.code + '; path=/;')
    else:
        response.set_header('Set-Cookie', 'anvioView=' + login + '|' + name + '; path=/;')

    return response


def set_project(request, userdb, response):
    set_default_headers(response)
    return json.dumps(userdb.set_project(get_user(request, userdb, response), request.forms.get('project')))


def get_current_project(request, userdb, response):
    set_default_headers(response)
    return json.dumps(userdb.get_current_project(request))


def get_current_project_files(request, userdb, response):
    set_default_headers(response)
    return json.dumps(userdb.get_current_project_files(request))


def get_current_project_archive(request, userdb, response):
    set_default_headers(response)

    zip_name, zip_archive = userdb.get_current_project_archive(request)

    # prepare response to download ZIP archive
    response.headers['Content-Type'] = 'application/zip, application/octet-stream; charset=UTF-8'
    response.headers['Content-Disposition'] = 'attachment; filename="{}"'.format(zip_name)

    return zip_archive.getvalue()


def debug(source, request):
    run.warning(None, header=source)
    run.warning(request.forms.dict, 'Forms', lc='yellow')
    run.warning(request.files.dict, 'Files', lc='yellow')
    run.warning(request.headers.items(), 'Headers', lc='yellow')
    run.warning(request.method, 'Method', lc='yellow')


def delete_project(request, userdb, response):
    set_default_headers(response)
    return json.dumps(userdb.delete_project(get_user(request, userdb, response), request.forms.get('project'), request.forms.get('user')))


def update_project(request, userdb, response):
    set_default_headers(response)
    return json.dumps(userdb.update_project(get_user(request, userdb, response), json.loads(request.body.readline())))


def share_project(request, userdb, response):
    set_default_headers(response)
    return json.dumps(userdb.create_view(get_user(request, userdb, response), request.forms.get('name'), request.forms.get('project'), request.forms.get('public')))


def delete_share(request, userdb, response):
    set_default_headers(response)
    return json.dumps(userdb.delete_view(get_user(request, userdb, response), request.forms.get('name')))


def receive_upload_file(request, userdb, response):
    set_default_headers(response)
    user = get_user(request, userdb, response)

    # check mandatory values
    if not user:
        return '{ "status": "error", "message": "you need to be logged in to create a project", "data": null }'

    if not request.files.get('treeFile'):
        return '{ "status": "error", "message": "you need to upload a tree file", "data": null }'

    if not request.forms.get('title'):
        return '{ "status": "error", "message": "a title is required to create a project", "data": null }'

    # create the project
    retval = userdb.create_project(user['data'], request.forms.get('title'), request.forms.get('description'))

    if not retval['status'] == 'ok':
        return json.dumps(retval)

    project = retval['data']

    # set user project to the new one
    retval = userdb.set_project(user['data'], project['name'])

    if not retval['status'] == 'ok':
        return json.dumps(retval)

    # save the uploaded files to the project directory
    basepath = userdb.users_data_dir + '/userdata/' + user['data']['path'] + '/' + project['path'] + '/'


    # tree and data fiels are mandatory
    request.files.get('treeFile').save(basepath + 'treeFile')

    if request.files.get('dataFile'):
        request.files.get('dataFile').save(basepath + 'dataFile')

    if request.files.get('fastaFile'):
        request.files.get('fastaFile').save(basepath + 'fastaFile')

    # check if we have samples information
    createSamplesDB = False
    samplesOrderPath = None
    if request.files.get('samplesOrderFile'):
        request.files.get('samplesOrderFile').save(basepath + 'samplesOrderFile')
        samplesOrderPath = basepath + 'samplesOrderFile'
        createSamplesDB = True

    samplesInfoPath = None
    if request.files.get('samplesInformationFile'):
        request.files.get('samplesInformationFile').save(basepath + 'samplesInformationFile')
        samplesInfoPath = basepath + 'samplesInformationFile'
        createSamplesDB = True

    # create a samples database if needed
    if createSamplesDB:
        sample = dbops.SamplesInformationDatabase(basepath + 'samples.db')
        try:
            sample.create(samplesInfoPath, samplesOrderPath)
        except Exception as e:
            return json.dumps({'status': 'error', 'message': "That one did not go as expected. Here is the error: %s" % e, "data": None})

    # all files are uploaded, do a sanity check
    retval = userdb.get_the_interactive_object(basepath, read_only=False)
    if not retval['status'] == 'ok':
        # if the files are not ok, the project needs to be deleted
        userdb.delete_project(user['data'], project['name'])

        # return the error to the ui
        return json.dumps(retval)

    return '{ "status": "ok", "message": "project created", "data": null }'


def receive_additional_upload_file(request, userdb, response):
    set_default_headers(response)
    if not request.files.get('uploadFile'):
        return '{ "status": "error", "message": "you did not upload a file", "data": null }'

    if not request.forms.get('project'):
        return '{ "status": "error", "message": "you did not specify a project", "data": null }'

    user = get_user(request, userdb, response)
    if not user:
        return '{ "status": "error", "message": "you need to be logged in to upload additional data", "data": null }'

    user = user['data']

    project = userdb.get_project(user, request.forms.get('project'))

    if not project['status'] == 'ok':
        return json.dumps(project)

    basepath = userdb.users_data_dir + '/userdata/' + user['path'] + '/' + project['data']['path'] + '/'

    fileType = 'additionalFile'
    validFileType = {'additionalFile': True, 'dataFile': True, 'treeFile': True, 'fastaFile': True, 'samplesOrderFile': True, 'samplesInformationFile': True}
    if request.forms.get('type'):
        fileType = request.forms.get('type')
        if not validFileType[fileType]:
            return '{ "status": "error", "message": "Invalid file type selected for upload", "data": null }'

    message = 'added'
    filePath = basepath + fileType
    if os.path.isfile(filePath):
        os.remove(filePath)
        message = 'updated'

    request.files.get('uploadFile').save(filePath)

    # check if this is a samples order or samples info file in which case the samples.db needs to be
    # created / updated
    if fileType == 'samplesOrderFile' or fileType == 'samplesInformationFile':
        # if there was a previous samples.db, it needs to be removed
        samplesDBPath = basepath + 'samples.db'
        if os.path.isfile(samplesDBPath):
            os.remove(samplesDBPath)

        # create a new samples.db
        samplesInfoPath = None
        samplesOrderPath = None
        if os.path.isfile(basepath + 'samplesInformationFile'):
            samplesInfoPath = basepath + 'samplesInformationFile'
        if os.path.isfile(basepath + 'samplesOrderFile'):
            samplesOrderPath = basepath + 'samplesOrderFile'
        if samplesInfoPath and samplesOrderPath:
            sample = dbops.SamplesInformationDatabase(samplesDBPath)
            try:
                sample.create(samplesInfoPath, samplesOrderPath)
            except Exception as e:
                return json.dumps({'status': 'error', 'message': "That one did not go as expected. Here is the error: %s" % e, "data": None})

    return '{ "status": "ok", "message": "file ' + message + '", "data": null }'

def admin_data(request, userdb, response):
    set_default_headers(response)
    user = get_user(request, userdb, response)
    if user and user['status'] == 'ok':
        if user['data']['clearance'] == 'admin':
            filterhash = {}
            fields = ['firstname', 'lastname', 'login', 'email', 'login', 'accepted', 'affiliation', 'clearance', 'date', 'visit']
            for field in fields:
                if field in request.query:
                    filterhash[field] = request.query[field]

            offset = 0
            limit = 25
            order = 'lastname'
            direction = 'ASC'
            if 'offset' in request.query:
                offset = request.query['offset']
            if 'limit' in request.query:
                limit = request.query['limit']
            if 'order' in request.query:
                order = request.query['order']
            if 'direction' in request.query:
                direction = request.query['direction']

            return json.dumps(userdb.user_list(offset, limit, order, direction, filterhash))
        else:
            return '{ "status": "error", "message": "You need to be an administrator to view this data", "data": null }'
    else:
        return '{ "status": "error", "message": "You need to be logged in to view this data", "data": null }'


def admin_project_data(request, userdb, response):
    set_default_headers(response)
    user = get_user(request, userdb, response)
    if user['status'] == 'ok':
        if user['data']['clearance'] == 'admin':
            filterhash = {}
            fields = ['name', 'user', 'path', 'description', 'views', 'metadata']
            for field in fields:
                if field in request.query:
                    filterhash[field] = request.query[field]

            offset = 0
            limit = 25
            order = 'name'
            direction = 'ASC'
            if 'offset' in request.query:
                offset = request.query['offset']
            if 'limit' in request.query:
                limit = request.query['limit']
            if 'order' in request.query:
                order = request.query['order']
            if 'direction' in request.query:
                direction = request.query['direction']

            return json.dumps(userdb.project_list(offset, limit, order, direction, filterhash))
        else:
            return '{ "status": "error", "message": "You need to be an administrator to view this data", "data": null }'
    else:
        return json.dumps(user)


def admin_project_details(request, userdb, response):
    set_default_headers(response)
    user = get_user(request, userdb, response)
    if user['status'] == 'ok':
        if user['data']['clearance'] == 'admin':
            return json.dumps(userdb.project_admin_details(request.query['project'], request.query['user']))
        else:
            return '{ "status": "error", "message": "You need to be an administrator to view this data", "data": null }'
    else:
        return json.dumps(user)
