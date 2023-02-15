import importlib.resources as pkg_resources
from . import yeast, human, mouse

__all__ = ['Resources']


class Resources:
    """An object to provide access to package resources to the user"""

    _configured_organisms = ['yeast', 'human', 'mouse']

    def __init__(self) -> None:
        self._yeast_resources = {
            'chr_map':
                pkg_resources.read_text(yeast, "chr_map.csv"),
            'barcode_details':
                pkg_resources.read_text(yeast, "barcode_details.json"),
            'background_sir4':
                pkg_resources.read_text(yeast, "background_dsir4.csv"),
            'background_adh1':
                pkg_resources.read_text(yeast, "background_adh1.csv"),
            'regions_not_orf':
                pkg_resources.read_text(yeast, "regions_not_orf.csv"),  # noqa
            'yiming_regions':
                pkg_resources.read_text(yeast, "yiming_regions.csv")  # noqa
        }
        self._human_resources = {
            'chr_map':
                pkg_resources.read_text(human, "chr_map.csv"),
            'barcode_details':
                pkg_resources.read_text(human, "barcode_details.json"),
        }
        self._mouse_resources = {
            'barcode_details':
                pkg_resources.read_text(mouse, "barcode_details.json"),
        }

    @property
    def configured_organisms(self):
        """list of organisms for which there are resources"""
        return self._configured_organisms

    @property
    def yeast_chr_map(self):
        """yeast chromosome name map"""
        return self._yeast_resources['chr_map']

    @property
    def yeast_background_sir4(self):
        """background hops from sir4 for yeast data in qBed format"""
        return self._yeast_resources['background_sir4']

    @property
    def yeast_background_adh1(self):
        """background hops from minus adh1 for yeast data in qBed format"""
        return self._yeast_resources['background_adh1']

    @property
    def yeast_promoters_yiming(self):
        """yeast promoters as defined by yiming"""
        return self._yeast_resources['yiming_regions']

    @property
    def yeast_promoters_not_orf(self):
        """Rob's original notOrf from the sacCer3 genome"""
        return self._yeast_resources['regions_not_orf']

    @property
    def yeast_not_orf(self):
        """Bed file of not ORF regions in sacCer3. 

        columns: Name, Chr, Start, Stop, Up, Down, size. 
        Up and Down are upstream/downstream features

        """
        return self._yeast_resources['notOrf_sacCer3.txt']

    @property
    def yeast_barcode_details(self):
        """An example barcode details file"""
        return self._yeast_resources['barcode_details']

    @property
    def human_resources(self):
        """Resources for human calling cards data"""
        return self._human_resources

    # def write_barcode_details(self,organism:str, dirpath:str) -> None:
    # 	"""Write a barcode details json to file

    # 	Args:
    # 		organism (str): organism for which to create an example barcode details file
    # 		dirpath (str): path to a directory to which to output the json
    # 	"""

    # 	if organism not in self.configured_organisms:
    # 		raise AttributeError(f"Recognized organisms are "\
    # 			f"{self.configured_organisms}. {organism} is not valid.")

    # 	if not os.path.isdir(dirpath):
    # 		raise NotADirectoryError(f"{dirpath} is not a valid directory.")

    #     with open(os.path.join(dirpath, f"{organism}_barcode_details_example.json", 'w') as f:
    #         f.write(cc_resources.yeast_barcode_details)


# @property
    # def mouse_resources(self):
    # 	"""Resources for human calling cards data"""
    # 	return self._human_resources
