import os
import ssl
import gzip
import socket
import shutil
import webbrowser
import urllib.request, urllib.error, urllib.parse

import Bio.PDB as PDB

import anvio
import anvio.constants as constants
import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError
from anvio.terminal import Run, Progress, SuppressAllOutput

from anvio.utils.algorithms import human_readable_file_size



def get_port_num(port_num = 0, ip='0.0.0.0', run=Run()):
    """Get a port number for the `ip` address."""

    try:
        port_num = int(port_num) if port_num else 0
    except Exception as e:
        raise ConfigError("Not a happy port number :/ %s." % e)

    if not port_num:
        port_num = get_next_available_port_num(constants.default_port_number)

        if not port_num:
            raise ConfigError("Anvi'o searched a bunch of port numbers starting from %d, but failed "
                               "to find an available one for you. Maybe you should specify one :/")
    else:
        if is_port_in_use(port_num):
            raise ConfigError("The port number %d seems to be in use :/" % port_num)

    if os.getuid() and port_num < 1024:
        run.warning("Using the port number %d requires superuser priviliges, which your user does not "
                    "seem to have. Since anvi'o does not know anything about your system configuraiton, "
                    "you are free to go for now. But be prepared for a failed attempt to use this port "
                    "number to serve stuff." % port_num)

    return port_num



def get_next_available_port_num(start=constants.default_port_number, look_upto_next_num_ports=100, ip='0.0.0.0'):
    """Starts from 'start' and incrementally looks for an available port
       until 'start + look_upto_next_num_ports', and returns the first
       available one."""
    for p in range(start, start + look_upto_next_num_ports):
        if not is_port_in_use(p, ip):
            return p

    return None



def is_port_in_use(port, ip='0.0.0.0'):
    in_use = False
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    result = sock.connect_ex((ip, port))

    if result == 0:
        in_use = True

    sock.close()
    return in_use



def download_file(url, output_file_path, check_certificate=True, progress=Progress(), run=Run()):
    filesnpaths.is_output_file_writable(output_file_path)

    if anvio.DEBUG:
        run.warning(None, header="DOWNLOADING FILE", overwrite_verbose=True, nl_before=1)
        run.info('Source URL', url, overwrite_verbose=True)
        run.info('Output path', output_file_path, overwrite_verbose=True, nl_after=1)

    try:
        if check_certificate:
            response = urllib.request.urlopen(url)
        else:
            response = urllib.request.urlopen(url, context=ssl._create_unverified_context())
    except Exception as e:
        raise ConfigError(f"Something went wrong with your download attempt. Here is the "
                          f"problem for the url {url}: '{e}'")

    file_size = 0
    if 'Content-Length' in response.headers:
        file_size = int(response.headers['Content-Length'])

    f = open(output_file_path, 'wb')

    progress.new('Downloading "%s"' % os.path.basename(output_file_path))
    progress.update('...')

    downloaded_size = 0
    counter = 0
    while True:
        buffer = response.read(10000)

        if buffer:
            downloaded_size += len(buffer)
            f.write(buffer)

            if counter % 500 == 0:
                if file_size:
                    progress.update('%.1f%%' % (downloaded_size * 100.0 / file_size))
                else:
                    progress.update('%s' % human_readable_file_size(downloaded_size))
        else:
            break

        counter += 1

    f.close()

    progress.end()
    run.info('Downloaded successfully', output_file_path)



def get_remote_file_content(url, gzipped=False, timeout=None):
    import requests
    from io import BytesIO

    if timeout:
        remote_file = requests.get(url, timeout=timeout)
    else:
        remote_file = requests.get(url)

    if remote_file.status_code == 404:
        raise ConfigError("Bad news. The remote file at '%s' was not found :(" % url)

    if gzipped:
        buf = BytesIO(remote_file.content)
        fg = gzip.GzipFile(fileobj=buf)
        return fg.read().decode('utf-8')

    return remote_file.content.decode('utf-8')



def get_anvio_news():
    """Reads news from anvi'o repository.

    The format of the news file is expected to be like this:

        # Title with spaces (01.01.1970) #
        Lorem ipsum, dolor sit amet
        ***
        # Title with spaces (01.01.1970) #
        Lorem ipsum, dolor sit amet
        ***
        # Title with spaces (01.01.1970) #
        Lorem ipsum, dolor sit amet

    Returns
    =======
    news : list
        A list of dictionaries per news item
    """

    try:
        news = get_remote_file_content(constants.anvio_news_url, timeout=1)
    except Exception as e:
        raise ConfigError(f"Something went wrong reading the anvi'o news :/ This is what the "
                          f"downstream library had to say: {e}")

    news_items = []
    for news_item in news.split('***'):
        if len(news_item) < 5:
            # too short to parse, just skip it
            continue

        news_items.append({'date': news_item.split("(")[1].split(")")[0].strip(),
                           'title': news_item.split("#")[1].split("(")[0].strip(),
                           'content': news_item.split("#\n")[1].strip()})

    return news_items



def download_protein_structure(protein_code, output_path=None, chain=None, raise_if_fail=True):
    """Downloads protein structures using Biopython.

    Parameters
    ==========
    protein_code : str
        Each element is a 4-letter protein code

    output_path : str
        Path where structure is written to. Temporary directory is chosen if None

    chain : str, None
        If None, all chains remain in the PDB file. If specified, only the chain with the chain ID
        `chain` will be saved.

    raise_if_fail : bool, True
        If the file does not download, raise an error

    Returns
    =======
    output : output_path
        Returns the filepath of the written file. Returns None if download failed
    """

    output_dir = os.path.dirname(output_path)
    if output_dir == '': output_dir = '.'

    pdb_list = PDB.PDBList()

    # NOTE This path is determined by Biopython's fn `pdb_list.retive_pdb_file`. If the logic in
    #      that function that determines the path name is changed, `download_protein_structure` will
    #      break because `temp_output_path` will be wrong.
    temp_output_path = os.path.join(output_dir, f"pdb{protein_code.lower()}.ent")

    try:
        with SuppressAllOutput():
            # We suppress output that looks like this:
            # >>> WARNING: The default download format has changed from PDB to PDBx/mmCif
            # >>> Downloading PDB structure '5w6y'...
            pdb_list.retrieve_pdb_file(protein_code, file_format='pdb', pdir=output_dir, overwrite=True)
    except:
        pass

    if not filesnpaths.is_file_exists(temp_output_path, dont_raise=True):
        # The file wasn't downloaded
        if raise_if_fail:
            raise ConfigError("The protein %s could not be downloaded. Are you connected to internet?" % protein_code)
        else:
            return None

    if chain is not None:
        class ChainSelect(PDB.Select):
            def accept_chain(self, chain_obj):
                return 1 if chain_obj.get_id() == chain else 0

        p = PDB.PDBParser()
        try:
            structure = p.get_structure(None, temp_output_path)
        except:
            # FIXME Something very rare happened on Biopython's end. We silently return the whole
            # file instead of only the chain. Here is one such reason for failure we stumbled upon:
            # https://github.com/biopython/biopython/issues/2819
            shutil.move(temp_output_path, output_path)
            return output_path

        # Overwrite file with chain-only structure
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(temp_output_path, ChainSelect())

    shutil.move(temp_output_path, output_path)

    return output_path



def open_url_in_browser(url, browser_path=None, run=Run()):
    if browser_path:
        filesnpaths.is_file_exists(browser_path)
        run.info_single('You are launching an alternative browser. Keep an eye on things!', mc='red', nl_before=1)
        webbrowser.register('users_preferred_browser', None, webbrowser.BackgroundBrowser(browser_path))
        webbrowser.get('users_preferred_browser').open_new(url)
    else:
        webbrowser.open_new(url)

