import time

import anvio.filesnpaths as filesnpaths

from anvio.errors import ConfigError, FilesNPathsError
from anvio.terminal import Run

from anvio.utils.misc import HTMLColorToRGB


def gen_gexf_network_file(units, samples_dict, output_file, sample_mapping_dict=None,
                               unit_mapping_dict=None, project=None, sample_size=8, unit_size=2,
                               skip_sample_labels=False, skip_unit_labels=False):
    """A function that generates an XML network description file for Gephi.

       Two minimum required inputs are `units`, and `samples_dict`.

       Simply, `samples_dict` is a dictionary that shows the distribution of `units` and their
       frequencies across samples. Here is an example `units` variable (which is a type of `list`):

            units = ['unit_1', 'unit_2', ... 'unit_n']

       and a corresponding `samples_dict` would look like this:

            samples_dict = {'sample_1': {'unit_1': 0.5,
                                        'unit_2': 0.2,
                                         ...,
                                         'unit_n': 0.1
                                        },
                            'sample_2': { (...)
                                            },
                            (...),
                            'sample_n': { (...)
                                            }
                        }
    """

    filesnpaths.is_output_file_writable(output_file)

    def RepresentsFloat(s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    output = open(output_file, 'w')

    samples = sorted(samples_dict.keys())
    sample_mapping_categories = sorted([k for k in list(sample_mapping_dict.values())[0].keys() if k != 'colors']) if sample_mapping_dict else None
    unit_mapping_categories = sorted([k for k in list(unit_mapping_dict.keys()) if k not in ['colors', 'labels']]) if unit_mapping_dict else None

    sample_mapping_category_types = []
    if sample_mapping_dict:
        for category in sample_mapping_categories:
            if RepresentsFloat(list(sample_mapping_dict.values())[0][category]):
                sample_mapping_category_types.append('double')
            else:
                sample_mapping_category_types.append('string')

    output.write('''<?xml version="1.0" encoding="UTF-8"?>\n''')
    output.write('''<gexf xmlns:viz="http:///www.gexf.net/1.1draft/viz" xmlns="http://www.gexf.net/1.2draft" version="1.2">\n''')
    output.write('''<meta lastmodifieddate="2010-01-01+23:42">\n''')
    output.write('''    <creator>Oligotyping pipeline</creator>\n''')
    if project:
        output.write('''    <creator>Network description for %s</creator>\n''' % (project))
    output.write('''</meta>\n''')
    output.write('''<graph type="static" defaultedgetype="undirected">\n\n''')

    if sample_mapping_dict:
        output.write('''<attributes class="node" type="static">\n''')
        for i in range(0, len(sample_mapping_categories)):
            category = sample_mapping_categories[i]
            category_type = sample_mapping_category_types[i]
            output.write('''    <attribute id="%d" title="%s" type="%s" />\n''' % (i, category, category_type))
        output.write('''</attributes>\n\n''')

    # FIXME: IDK what the hell is this one about:
    if unit_mapping_dict:
        output.write('''<attributes class="edge">\n''')
        for i in range(0, len(unit_mapping_categories)):
            category = unit_mapping_categories[i]
            output.write('''    <attribute id="%d" title="%s" type="string" />\n''' % (i, category))
        output.write('''</attributes>\n\n''')

    output.write('''<nodes>\n''')
    for sample in samples:
        if skip_sample_labels:
            output.write('''    <node id="%s">\n''' % (sample))
        else:
            output.write('''    <node id="%s" label="%s">\n''' % (sample, sample))

        output.write('''        <viz:size value="%d"/>\n''' % sample_size)

        if sample_mapping_dict and 'colors' in sample_mapping_dict[sample]:
            output.write('''        <viz:color r="%d" g="%d" b="%d" a="1"/>\n''' %\
                                             HTMLColorToRGB(sample_mapping_dict[sample]['colors'], scaled=False))

        if sample_mapping_categories:
            output.write('''        <attvalues>\n''')
            for i in range(0, len(sample_mapping_categories)):
                category = sample_mapping_categories[i]
                output.write('''            <attvalue id="%d" value="%s"/>\n''' % (i, sample_mapping_dict[sample][category]))
            output.write('''        </attvalues>\n''')

        output.write('''    </node>\n''')

    for unit in units:
        if skip_unit_labels:
            output.write('''    <node id="%s">\n''' % (unit))
        else:
            if unit_mapping_dict and 'labels' in unit_mapping_dict:
                output.write('''    <node id="%s" label="%s">\n''' % (unit, unit_mapping_dict['labels'][unit]))
            else:
                output.write('''    <node id="%s">\n''' % (unit))
        output.write('''        <viz:size value="%d" />\n''' % unit_size)

        if unit_mapping_categories:
            output.write('''        <attvalues>\n''')
            for i in range(0, len(unit_mapping_categories)):
                category = unit_mapping_categories[i]
                output.write('''            <attvalue id="%d" value="%s"/>\n''' % (i, unit_mapping_dict[category][unit]))
            output.write('''        </attvalues>\n''')

        output.write('''    </node>\n''')

    output.write('''</nodes>\n''')

    edge_id = 0
    output.write('''<edges>\n''')
    for sample in samples:
        for i in range(0, len(units)):
            unit = units[i]
            if samples_dict[sample][unit] > 0.0:
                if unit_mapping_dict:
                    output.write('''    <edge id="%d" source="%s" target="%s" weight="%f">\n''' % (edge_id, unit, sample, samples_dict[sample][unit]))
                    if unit_mapping_categories:
                        output.write('''        <attvalues>\n''')
                        for i in range(0, len(unit_mapping_categories)):
                            category = unit_mapping_categories[i]
                            output.write('''            <attvalue id="%d" value="%s"/>\n''' % (i, unit_mapping_dict[category][unit]))
                        output.write('''        </attvalues>\n''')
                    output.write('''    </edge>\n''')
                else:
                    output.write('''    <edge id="%d" source="%s" target="%s" weight="%f" />\n''' % (edge_id, unit, sample, samples_dict[sample][unit]))


                edge_id += 1
    output.write('''</edges>\n''')
    output.write('''</graph>\n''')
    output.write('''</gexf>\n''')

    output.close()



def run_selenium_and_export_svg(url, output_file_path, browser_path=None, run=Run()):
    if filesnpaths.is_file_exists(output_file_path, dont_raise=True):
        raise FilesNPathsError("The output file already exists. Anvi'o does not like overwriting stuff.")

    filesnpaths.is_output_file_writable(output_file_path)

    try:
        from selenium import webdriver
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait
        from selenium.webdriver.support import expected_conditions as EC
        from selenium.common.exceptions import TimeoutException
    except:
        raise ConfigError("You want to export SVGs? Well, you need the Python library 'selenium' to be able to "
                          "do that but you don't have it. If you are lucky, you probably can install it by "
                          "typing 'pip install selenium' or something :/")

    if browser_path:
        filesnpaths.is_file_exists(browser_path)
        run.info_single('You are launching an alternative browser. Keep an eye on things!', mc='red', nl_before=1)
        driver = webdriver.Chrome(executable_path=browser_path)
    else:
        driver = webdriver.Chrome()

    driver.wait = WebDriverWait(driver, 10)
    driver.set_window_size(1920, 1080)
    driver.get(url)

    try:
        WebDriverWait(driver, 300).until(EC.text_to_be_present_in_element((By.ID, "title-panel-second-line"), "Current view"))
    except TimeoutException:
        print("Timeout occured, could not get the SVG drawing in 600 seconds.")
        driver.quit()
    time.sleep(1)

    driver.execute_script("exportSvg(true);")
    time.sleep(1)

    svg = driver.find_element_by_id('panel-center')

    svg_file = open(output_file_path, 'w')
    svg_file.write(svg.get_attribute('innerHTML'))
    svg_file.close()
    driver.quit()

    run.info_single('\'%s\' saved successfully.' % output_file_path)

