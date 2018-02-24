from configparser import ConfigParser
from dataPreprocessing.configValues import ConfigValues


def parseConfig(fileName):
    print("Parse config", fileName)
    config = ConfigParser()
    config.read(fileName)
    return ConfigValues(config)