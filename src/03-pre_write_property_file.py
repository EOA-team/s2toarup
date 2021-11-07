
"""

@author:    Lukas Graf (D-USYS, ETHZ)

@purpose:   This scripts generates xml parameter files
            to run the S2 L1C-RUT for the single uncertainty
            contributors.

"""

from pathlib import Path
from typing import Optional
import xml.etree.ElementTree as ET


def write_properties(
        out_dir: Path,
        param_xml_template: Optional[Path]='parameters.xml'
    ) -> None:
    """
    Creates the property files required to run L1C-RUT
    for each uncertainty contributor. I.e., there will be as
    many files as uncertainty contributors, each named
    according to the uncertainty contributor

    :param out_dir:
        directory where to save the xml files to
    :param param_xml_template:
        xml file serving as template for creating the single
        xml parameter files
    """

    # parse the template xml file first
    tree = ET.parse(param_xml_template)
    root = tree.getroot()

    # by default, all uncertainty contributors (i.e., children of the
    # root element) are set to False
    params = [child for child in root]

    # now we can create n xml files where n is the number of uncertainty
    # contributors (i.e., number of children)
    for idx, param in enumerate(params):

        # create output file
        fname = out_dir.joinpath(f'{param.tag}.properties')
        with open(fname, 'w+') as dst:
            dst.writelines(f'{param.tag}=true' + '\n')

            # leave the other parameters unchanged (set to False)
            for jdx, param2 in enumerate(params):
                if jdx == idx: continue
                dst.writelines(f'{param2.tag}=false' + '\n')

        # write xml tree to file
        print(f'Wrote property file for "{param.tag}" uncertainty contributor')


if __name__ == '__main__':
    
    out_dir = Path('./../S2A_MSIL1C_orig')
    write_properties(out_dir=out_dir)
