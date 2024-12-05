import os
import sys
from pathlib import Path
from multiprocessing import Pool, cpu_count, get_context
from functools import partial
import logging

from copy import deepcopy

# Configure Logging
logging.basicConfig(
    level=logging.INFO,  # Change to DEBUG for more detailed logs
    format='%(asctime)s [%(levelname)s] [Process %(process)d] %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),  # Logs to console
        logging.FileHandler("file_generator.log")  # Logs to a file
    ]
)

logger = logging.getLogger(__name__)

# Update sys.path for module imports
sys.path.extend([
    '../AEM',
    '../MEM',
    '../file-conversion',
    '../../jsatorb-visibility-service/src'
])

# Define the placeholder conversion function
def ccsds2cic_in_memory(ccsds_data, body):
    """
    Converts CCSDS data to CIC format in-memory without writing to disk.
    Placeholder function: Implement actual in-memory conversion logic.
    """
    # TODO: Implement in-memory CCSDS to CIC conversion
    # For demonstration, returning the original data
    return ccsds_data

class FileGenerator:
    """ Master class that calls MissionAnalysis.py, MEMGenerator.py, and AEMGenerator.py,
    to generate OEM, AEM, and MEM files for given satellites, time options, and central body.
    It also generates a color file for VTS"""

    def __init__(self, stringDateStart, stringDateEnd, step, stringBody, satellites, groundStations, options):
        self.step = step
        self.body = stringBody
        self.satellites = satellites
        self.groundStations = groundStations
        self.options = options

        # Store parameters for worker initialization
        self.stringDateStart = stringDateStart
        self.stringDateEnd = stringDateEnd

    def _initialize_worker(self):
        """
        Initializes Orekit within the worker process.
        """
        try:
            import orekit
            from orekit.pyhelpers import setup_orekit_curdir
            orekit.initVM()
            setup_orekit_curdir()
            #logger.info("JVM initialized in worker process.")
        except Exception as e:
            logger.error(f"Failed to initialize JVM in worker: {e}")
            raise

    def _process_satellite(self, sat, nameFolder):
        """
        Processes a single satellite: generates OEM, AEM, MEM, and Color files.
        Initializes Orekit within the worker process to avoid JVM conflicts.
        """
        try:
            logger.info(f"Starting processing of satellite: {sat['name']}")

            # Initialize JVM in worker
            self._initialize_worker()

            # Import Java classes after initializing Orekit
            from MissionAnalysis import HAL_MissionAnalysis
            from AEMGenerator import AEMGenerator
            from MEMGenerator import MEMGenerator
            from ColorGenerator import ColorGenerator
            from ccsds2cic import ccsds2cic

            local_options = deepcopy(self.options)
            responses = {}

            # OEM Generation (CARTESIAN)
            if 'CARTESIAN' in local_options:
                #logger.info(f"Generating OEM for satellite: {sat['name']}")
                newMission = HAL_MissionAnalysis(self.step, self.stringDateEnd, self.body)
                newMission.setStartTime(self.stringDateStart)
                newMission.addSatellite(sat)
                newMission.propagate()
                oem_data = newMission.getOEMEphemerids()
                # Convert OEM data from CCSDS to CIC in-memory
                oem_converted = ccsds2cic_in_memory(oem_data, self.body)
                responses['OEM'] = oem_converted
                oem_file = f"{nameFolder}/{sat['name']}_OEM_POSITION.TXT"
                oem_file_ccsds = f"{nameFolder}/{sat['name']}_OEM_POSITION.TXT_ccsds"
                with open(oem_file_ccsds, 'w') as file:
                    file.write(oem_data)
                ccsds2cic(oem_file_ccsds, oem_file, self.body)
                # Remove 'CARTESIAN' from local_options to avoid re-processing
                del local_options['CARTESIAN']
                os.remove(oem_file_ccsds)
                #logger.info(f"OEM generation completed for satellite: {sat['name']}")

            # AEM Generation (ATTITUDE)
            if 'ATTITUDE' in local_options:
                #logger.info(f"Generating AEM for satellite: {sat['name']}")
                aemGenerator = AEMGenerator(self.stringDateStart, self.step, self.stringDateEnd, self.body)
                aem_file = f"{nameFolder}/{sat['name']}_AEM_ATTITUDE.TXT"
                aem_file_ccsds = f"{nameFolder}/{sat['name']}_AEM_ATTITUDE.TXT_ccsds"
                aemGenerator.setSatellite(sat)
                aemGenerator.setAttitudeLaw(local_options['ATTITUDE'])
                aemGenerator.setFile(aem_file_ccsds)
                aem_data = aemGenerator.propagate()
                # Convert AEM data from CCSDS to CIC in-memory
                aem_converted = ccsds2cic_in_memory(aem_data, self.body)
                ccsds2cic(aem_file_ccsds, aem_file, self.body)
                responses['AEM'] = aem_converted
                os.remove(aem_file_ccsds)
                # Remove 'ATTITUDE' from local_options to avoid re-processing
                del local_options['ATTITUDE']
                #logger.info(f"AEM generation completed for satellite: {sat['name']}")

            # MEM Generation
            #logger.info(f"Generating MEM for satellite: {sat['name']}")
            memGenerator = MEMGenerator(self.stringDateStart, self.step, self.stringDateEnd, self.body)
            memGenerator.setSatellite(sat)
            for memType, memOption in local_options.items():
                if memType == 'VISIBILITY':
                    for gs in self.groundStations:
                        mem_file = f"{nameFolder}/{sat['name']}_MEM_VISIBILITY_{gs['name']}.TXT"
                        memGenerator.addMemVisibility(mem_file, gs)
                        #logger.info(f"Added MEM Visibility for ground station {gs['name']} to satellite {sat['name']}")
                else:
                    mem_file = f"{nameFolder}/{sat['name']}_MEM_{memType}.TXT"
                    memGenerator.addMemType(memType, mem_file)
                    #logger.info(f"Added MEM type {memType} to satellite {sat['name']}")
            memGenerator.propagate()
            responses['MEM'] = "MEM propagation completed."
            #logger.info(f"MEM generation completed for satellite: {sat['name']}")

            # Color Generation
            #logger.info(f"Generating Color file for satellite: {sat['name']}")
            color_file = f"{nameFolder}/{sat['name']}_COLOR.TXT"
            colorGenerator = ColorGenerator(self.stringDateStart, color_file, sat)
            colorGenerator.generate()
            #logger.info(f"Color file generated for satellite: {sat['name']}")

            if 'VISIBILITY' in self.options:
                listVisibilities = [
                    f"{nameFolder}/{sat['name']}_MEM_VISIBILITY_{gs['name']}.TXT"
                    for gs in self.groundStations
                ]
                colorGenerator.addVisibilities(listVisibilities)
                #logger.info(f"Added visibilities to Color file for satellite: {sat['name']}")

            logger.info(f"Completed processing of satellite: {sat['name']}")
            return responses

        except Exception as e:
            logger.error(f"Error processing satellite {sat['name']}: {e}")
            return {"Error": str(e)}

    def generate(self, nameFolder):
        """
        Generates files for all satellites using multiprocessing.
        """
        logger.info("Starting file generation process.")

        # Ensure necessary directories exist
        Path(nameFolder).mkdir(parents=True, exist_ok=True)
        data_folder = Path(nameFolder) / 'Data'
        data_folder.mkdir(parents=True, exist_ok=True)

        # Define a partial function for multiprocessing
        func = partial(self._process_satellite, nameFolder=nameFolder)

        # Use 'spawn' start method to avoid inheriting the main process's state
        ctx = get_context("spawn")

        try:
            with ctx.Pool(processes=cpu_count()) as pool:
                results = pool.map(func, self.satellites)
            logger.info("File generation completed successfully.")
            return results
        except KeyboardInterrupt:
            logger.warning("KeyboardInterrupt received. Terminating pool.")
            pool.terminate()
            pool.join()
            sys.exit(1)
        except Exception as e:
            logger.error(f"Multiprocessing pool encountered an error: {e}")
            pool.terminate()
            pool.join()
            raise

# Example Usage
if __name__ == "__main__":
    # Define your options, general settings, ground stations, and satellites
    options = {
        "CARTESIAN": {},
        "KEPLERIAN": {},
        "ATTITUDE": {
            "law": "LOF_LVLH"
        },
        "ECLIPSE": {}
    }

    general = {
        'timeStart': "2011-12-01T16:43:45",
        'timeEnd': "2011-12-02T16:43:45",
        'celestialBody': 'EARTH'
    }

    groundStations = [
        {
            "name": "isae",
            "latitude": 43,
            "longitude": 1.5,
            "altitude": 150,
            "elevation": 12
        },
        {
            "name": "cayenne",
            "latitude": 4.5,
            "longitude": -52.9,
            "altitude": 0,
            "elevation": 12
        }
    ]

    satellites = [
        {
            "name": "KepSat",
            "type": "keplerian",
            "sma": 7000000,
            "ecc": 0.007014455530245822,
            "inc": 51,
            "pa": 0,
            "raan": 0,
            "meanAnomaly": 0
        },
        {
            "name": "CartSat",
            "type": "cartesian",
            "x": -6142438.668,
            "y": 3492467.560,
            "z": -25767.25680,
            "vx": 505.8479685,
            "vy": 942.7809215,
            "vz": 7435.922231
        }
    ]

    # Path to new folder
    nameFolder = 'generated_files/'

    generator = FileGenerator(
        stringDateStart=general['timeStart'],
        stringDateEnd=general['timeEnd'],
        step=60,
        stringBody=general['celestialBody'],
        satellites=satellites,
        groundStations=groundStations,
        options=options
    )
    try:
        results = generator.generate(nameFolder)
        logger.info(f"Generation Results: {results}")
    except Exception as e:
        logger.error(f"File generation failed: {e}")
