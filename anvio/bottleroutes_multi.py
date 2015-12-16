# -*- coding: utf-8
"""
    Multi-User routes for bottle web server.

    Functions defined here are used from wihtin programs such as
    anvi-interactive, or anvi-refine.
"""

import os
import json
import anvio
import anvio.dbops as dbops

from bottle import redirect

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

    
def get_user_by_token(request, userdb, response):
    set_default_headers(response)
    retval = userdb.get_user_for_token(request.forms.get('token'))
        
    if retval[0]:
        del retval[1]['password']
        del retval[1]['path']
        del retval[1]['accepted']
        
    return json.dumps(retval)


def request_account(request, userdb, response):
    set_default_headers(response)
    retval = userdb.create_user(request.forms.get('firstname'), request.forms.get('lastname'), request.forms.get('email'), request.forms.get('login'), request.forms.get('password'))

    set_default_headers(response)
    if retval[0]:
        return '{ "OK": "'+retval[1]+'" }'
    else:
        return '{ "ERROR": "'+retval[1]+'" }'


def accept_user(request, userdb, response):
    retval = userdb.accept_user(request.query.login, request.query.code)
    if retval[0]:
        redirect('/app/accountOK.html')
    else:
        redirect('/app/accountBAD.html')
    
    return retval[1]


def reset_password(request, userdb, response):
    set_default_headers(response)
    if request.forms.get('email'):
        user = userdb.get_user_for_email(request.forms.get('email'))
        if user:
            userdb.reset_password(user)
            return '{ "OK": "password reset" }'
        else:
            return '{ "ERROR": "email address not found" }'
    else:
        return '{ "ERROR":"you need to pass an email address" }'


def check_availability(request, userdb, response):
    set_default_headers(response)
    if request.forms.get('email'):
        user = userdb.get_user_for_email(request.forms.get('email'))
        if user:
            return '{ "email": "this email is already taken" }'
        else:
            return '{ "email": "ok" }'
    elif request.forms.get('login'):
        user = userdb.get_user_for_login(request.forms.get('login'))
        if user:
            return '{ "login": "this login is already taken" }'
        else:
            return '{ "login": "ok" }'

        
def change_password(request, userdb, response):
    set_default_headers(response)
    if request.forms.get('login'):
        user = userdb.get_user_for_login(request.forms.get('login'))
        if user:
            if request.forms.get('password'):
                userdb.change_password(user, request.forms.get('password'))
                return '{ "OK": "password changed" }'
            else:
                return '{ "ERROR": "you must pass a password" }'
        else:
            return '{ "ERROR": "user not found" }'
    else:
        return '{ "ERROR": "you must provide a login" }'

    
def login_to_app(request, userdb, response):
    set_default_headers(response)
    retval = userdb.login_user(request.forms.get('login'), request.forms.get('password'))
    if retval[0]:
        response.set_header('Set-Cookie', 'anvioSession='+retval[1]["token"]+'; path=/; max-age='+str(60 * 60 * 24 * 14))
        
    return json.dumps(retval)


def logout_from_app(request, userdb, response):
    userdb.logout_user(request.forms.get('login'))
    return 'OK'


def set_view_cookie(request, userdb, response):
    if request.query.name:
        name = request.query.name
        token = request.query.code or ""
        response.set_header('Set-Cookie', 'anvioView='+name+'|'+token+'; path=/;')
        
        redirect('/app/index.html')
        

def set_project(request, userdb, response):
    set_default_headers(response)
    retval = get_user(request, userdb, response)
    if retval[0]:
        if request.forms.get('project'):
            userdb.set_project(retval[1]['login'], request.forms.get('project'))
            redirect('/app/index.html')
        else:
            return '{ "ERROR": "You need to specify a project name" }'
    else:
        return '{ "ERROR": "' + retval[1] + '" }'


def delete_project(request, userdb, response):
    set_default_headers(response)
    retval = get_user(request, userdb, response)
    if retval[0]:
        if request.forms.get('project'):
            userdb.delete_project(retval[1]['login'], request.forms.get('project'))
            redirect('/app/index.html')
        else:
            return '{ "ERROR": "You need to specify a project name" }'
    else:
        return '{ "ERROR": "' + retval[1] + '" }'
    

def share_project(request, userdb, response):
    set_default_headers(response)
    if not request.forms.get('name'):
        return '{ "ERROR": "no name specified for the share" }'
    
    if not request.forms.get('project'):
        return '{ "ERROR": "no project specified for the share" }'

    if not re.match("^[A-Za-z0-9_-]+$", request.forms.get('name')):
        return '{ "ERROR": "the share name contains invalid characters" }'
    
    retval = get_user(request, userdb, response)
    if not retval[0]:
        return '{ "ERROR": "no user logged in" }'

    if not retval[1]['project']:
        return '{ "ERROR": "no project selected" }'
        
    share = userdb.create_view(retval[1]['login'], request.forms.get('name'), request.forms.get('project'), request.forms.get('public'))
    if share[0]:
        return '{ "token": "'+share[1]+'", "project": "'+request.forms.get('project')+'", "name": "'+request.forms.get('name')+'", "public": "'+request.forms.get('public')+'"}'
    else:
        return '{ "ERROR": "'+share[1]+'" }'

    
def delete_share(request, userdb, response):
    set_default_headers(response)
    if not request.forms.get('name'):
        return '{ "ERROR": "no name specified for the share to delete" }'
    
    if not request.forms.get('project'):
        return '{ "ERROR": "no project specified for the share to delete" }'
    
    retval = get_user(request, userdb, response)
    if not retval[0]:
        return '{ "ERROR": "no user logged in" }'
        
    share = userdb.delete_view(retval[1]['login'], request.forms.get('name'))
    if share[0]:
        return '{ "OK": "view deleted", "project": "'+request.forms.get('project')+'", "view": "'+request.forms.get('project')+'" }'
    else:
        return '{ "ERROR": "'+share[1]+'" }'
    
    
def receive_upload_file(request, userdb, response):
    set_default_headers(response)

    retval = get_user(request, userdb, response)
    if not retval[0]:
        return '{ "ERROR": "you need to be logged in to create a project" }'
    
    if not request.forms.get('title'):
        return '{ "ERROR": "a title is required to create a project" }'

    if not request.files.get('treeFile'):
        return '{ "ERROR": "you need to upload a tree file" }'
    
    user = retval[1]
    retval = userdb.create_project(user['login'], request.forms.get('title'))

    if not retval[0]:
        return '{ "ERROR": "'+retval[1]+'" }'

    project = retval[1]

    retval = userdb.set_project(user['login'], project['name'])

    if not retval[0]:
        return '{ "ERROR": "'+retval[1]+'" }'
    
    basepath = userdb.users_data_dir + '/userdata/'+user['path']+'/'+project['path']+'/'
    
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
    profiledb = profile.create({"db_type": "profile", "contigs_db_hash": None});
        
    redirect('/app/index.html')


def receive_additional_upload_file(request, userdb, response):
    set_default_headers(response)
    if not request.files.get('additionalFile'):
        return '{ "ERROR": "you did not upload a file" }'

    if not request.forms.get('project'):
        return '{ "ERROR": "you did not specify a project" }'

    retval = get_user(request, userdb, response)
    if not retval[0]:
        return '{ "ERROR": "you need to be logged in to upload additional data" }'
    
    user = retval[1]

    if not retval[0]:
        return '{ "ERROR": "'+retval[1]+'" }'

    project = userdb.get_project(user['login'], request.forms.get('project'))

    basepath = userdb.users_data_dir + '/userdata/'+user['path']+'/'+project['path']+'/'
    
    request.files.get('additionalFile').save(basepath + 'additionalFile')
        
    return '{ "OK": "file added" }'

def admin_page(request, userdb, response):
    return True
