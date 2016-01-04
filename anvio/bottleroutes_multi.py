# -*- coding: utf-8
"""
    Multi-User routes for bottle web server.

    Functions defined here are used from wihtin programs such as
    anvi-interactive, or anvi-refine.
"""

import json
from bottle import redirect

import anvio.dbops as dbops


def set_default_headers(response):
    response.set_header('Content-Type', 'application/json')
    response.set_header('Pragma', 'no-cache')
    response.set_header('Cache-Control', 'no-cache, no-store, max-age=0, must-revalidate')
    response.set_header('Expires', 'Thu, 01 Dec 1994 16:00:00 GMT')
    

def get_user(request, userdb, response):
    # check if we have a cookie
    if request.get_cookie('anvioSession'):

        # we have a cookie, check if it is valid
        return userdb.get_user_for_token(request.get_cookie('anvioSession'))

    else:
        return False
    
def get_user_by_token(request, userdb, response):
    set_default_headers(response)
    return json.dumps(userdb.get_user_for_token(request.forms.get('token')))


def impersonate(request, userdb, response):
    set_default_headers(response)
    retval = get_user(request, userdb, response)
    if retval:
        if retval['data']['clearance'] == 'admin':
            imperson = userdb.get_user_for_login(request.forms.get('login'))
            if imperson['status'] == 'ok':
                response.set_header('Set-Cookie', 'anvioSession='+imperson['data']["token"]+'; path=/; max-age='+str(60 * 60 * 24 * 14))
            return json.dumps(imperson)


def request_account(request, userdb, response):
    set_default_headers(response)
    return json.dumps(userdb.create_user(request.forms.get('firstname'), request.forms.get('lastname'), request.forms.get('email'), request.forms.get('login'), request.forms.get('password'), request.forms.get('affiliation'), request.environ.get('REMOTE_ADDR'), ))


def accept_user(request, userdb, response):
    retval = userdb.accept_user(request.query.login, request.query.code)
    if retval['status'] == 'ok':
        redirect('/app/accountOK.html')
    else:
        redirect('/app/accountBAD.html')


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
        response.set_header('Set-Cookie', 'anvioSession='+retval['data']["token"]+'; path=/; max-age='+str(60 * 60 * 24 * 14))
        
    return json.dumps(retval)


def logout_from_app(request, userdb, response):
    return json.dumps(userdb.logout_user(request.forms.get('login')))


def set_view_cookie(request, userdb, response):
    if request.query.name:
        name = request.query.name
        token = request.query.code or ""
        response.set_header('Set-Cookie', 'anvioView='+name+'|'+token+'; path=/;')
        
        redirect('/app/index.html')
        

def set_project(request, userdb, response):
    set_default_headers(response)
    return json.dumps(userdb.set_project(get_user(request, userdb, response), request.forms.get('project')))


def delete_project(request, userdb, response):
    set_default_headers(response)
    return json.dumps(userdb.delete_project(get_user(request, userdb, response), request.forms.get('project')))

def share_project(request, userdb, response):
    set_default_headers(response)
    return json.dumps(share = userdb.create_view(get_user(request, userdb, response), request.forms.get('name'), request.forms.get('project'), request.forms.get('public')))

    
def delete_share(request, userdb, response):
    set_default_headers(response)
    return json.dumps(userdb.delete_view(get_user(request, userdb, response), request.forms.get('name')))
    
    
def receive_upload_file(request, userdb, response):
    set_default_headers(response)
    user = get_user(request, userdb, response)
    if not user:
        return '{ "status": "error", "message": "you need to be logged in to create a project", "data": null }'
    
    if not request.forms.get('title'):
        return '{ "status": "error", "message": "a title is required to create a project", "data": null }'

    if not request.files.get('treeFile'):
        return '{ "status": "error", "message": "you need to upload a tree file", "data": null }'
    
    retval = userdb.create_project(user['data'], request.forms.get('title'))

    if not retval['status'] == 'ok':
        return json.dumps(retval)

    project = retval['data']

    retval = userdb.set_project(user['data'], project['name'])

    if not retval['status'] == 'ok':
        return json.dumps(retval)
    
    basepath = userdb.users_data_dir + '/userdata/'+user['data']['path']+'/'+project['path']+'/'
    
    request.files.get('treeFile').save(basepath + 'treeFile')
    if request.files.get('fastaFile'):
        request.files.get('fastaFile').save(basepath + 'fastaFile')
    else:
        open(basepath + 'fastaFile', 'a').close()
    if request.files.get('dataFile'):
        request.files.get('dataFile').save(basepath + 'dataFile')
    else:
        open(basepath + 'dataFile', 'a').close()

    # create a profile database
    profile = dbops.ProfileDatabase(basepath + 'profile.db')
    profiledb = profile.create({"db_type": "profile", "contigs_db_hash": None})

    # create a samples database
    sample = dbops.SamplesInformationDatabase(basepath + 'samples.db')
    samplesdb = sample.create()
        
    redirect('/app/index.html')


def receive_additional_upload_file(request, userdb, response):
    set_default_headers(response)
    if not request.files.get('additionalFile'):
        return '{ "status": "error", "message": "you did not upload a file", "data": null }'

    if not request.forms.get('project'):
        return '{ "status": "error", "message": "you did not specify a project", "data": null }'

    user = get_user(request, userdb, response)
    if not user:
        return '{ "status": "error", "message": "you need to be logged in to upload additional data", "data": null }'
    
    user = user[data]

    project = userdb.get_project(user, request.forms.get('project'))

    basepath = userdb.users_data_dir + '/userdata/'+user['path']+'/'+project['path']+'/'
    
    request.files.get('additionalFile').save(basepath + 'additionalFile')
        
    return '{ "status": "ok", "message": "file added", "data": null }'

def admin_data(request, userdb, response):
    set_default_headers(response)
    user = get_user(request, userdb, response)
    if user['status'] == 'ok':
        if user['data']['clearance'] == 'admin':
            filterhash = {}
            fields = [ 'firstname', 'lastname', 'login', 'email', 'login', 'accepted', 'affiliation', 'clearance', 'date' ]
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
        return json.dumps(user)
