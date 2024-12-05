from distutils.dir_util import copy_tree
import os
import sys
import io
import zipfile
from pathlib import Path
from datetime import datetime
from bottle import Bottle, request, response, HTTPResponse

# Initialize the Bottle application
app = Bottle()

import json
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

#### Give visibility on processing modules called from the REST API
# Add mission analysis module
sys.path.append('../jsatorb-visibility-service/src')
# Add eclipses module
sys.path.append('../jsatorb-eclipse-service/src')
# Add Date conversion module
sys.path.append('../jsatorb-date-conversion/src')
# Add JSatOrb common module: AEM and MEM generators
sys.path.append('../jsatorb-common/src')
sys.path.append('../jsatorb-common/src/AEM')
sys.path.append('../jsatorb-common/src/MEM')
sys.path.append('../jsatorb-common/src/VTS')
# Add file conversion module
sys.path.append('../jsatorb-common/src/file-conversion')
# Add JSatOrb common module: Mission Data management
sys.path.append('../jsatorb-common/src/mission-mgmt')
# Add Constellation generator module
sys.path.append('../jsatorb-common/src/constellation')
# Add Coverage module
sys.path.append('../jsatorb-coverage-service/src')

import bottle
from bottle import request, response
from MissionAnalysis import HAL_MissionAnalysis
from WalkerConstellation import WalkerConstellation
from DateConversion import HAL_DateConversion
from EclipseCalculator import HAL_SatPos, EclipseCalculator
from FileGenerator import FileGenerator
from VTSGenerator import VTSGenerator
from ccsds2cic import ccsds2cic
from MissionDataManager import writeMissionDataFile, loadMissionDataFile, listMissionDataFile, duplicateMissionDataFile, \
    isMissionDataFileExists, deleteMissionDataFile
from CoverageGenerator import CoverageGenerator
from VTSGeneratorCoverage import VTSGeneratorCoverage
from datetime import datetime
import json


# Helper Functions

def enable_cors(fn):
    def _enable_cors(*args, **kwargs):
        response.headers['Access-Control-Allow-Origin'] = '*'
        response.headers['Access-Control-Allow-Methods'] = 'GET, POST, PUT, OPTIONS, DELETE'
        response.headers[
            'Access-Control-Allow-Headers'] = 'Origin, Accept, Content-Type, X-Requested-With, X-CSRF-Token'
        if request.method != 'OPTIONS':
            return fn(*args, **kwargs)

    return _enable_cors


def show_request(req):
    """Log received HTTP request."""
    logger.info("RECEIVED REQUEST --------------------------------------------------")
    logger.info(req)
    logger.info("END OF RECEIVED REQUEST -------------------------------------------")


def show_response(res):
    """Log sent HTTP response."""
    logger.info("SENT RESPONSE (truncated to 1000 char) ----------------------------")
    logger.info(res[:1000])
    logger.info("END OF SENT RESPONSE ----------------------------------------------")


def bool_to_rest_status(value):
    """Convert a boolean to a REST status value {"SUCCESS", "FAIL"}."""
    return "SUCCESS" if value else "FAIL"


def build_smd_response(status, message, data):
    """
    Build a formatted REST response as a dictionary:
    {
        "status": <operation status: "SUCCESS" or "FAIL">,
        "message": <error message if "FAIL" is returned>,
        "data": <response data>
    }
    """
    return {"status": status, "message": message, "data": data}


def error_response(error_name):
    """Return an error response as JSON."""
    return json.dumps({"error": error_name})


def get_list_of_files(dir_name):
    """
    Efficiently retrieve all file paths within a directory tree.
    """
    return [str(path) for path in Path(dir_name).rglob('*') if path.is_file()]


def zipped_vts_response(vts_folder, mission):
    """
    Create an HTTP response containing the VTS compressed data structure.
    """
    buf = io.BytesIO()
    list_of_files = get_list_of_files(vts_folder)

    with zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED) as zipfh:
        for individual_file in list_of_files:
            # Compute relative path for the archive
            arcname = os.path.relpath(individual_file, vts_folder)
            zipfh.write(individual_file, arcname)

    buf.seek(0)
    filename = f'vts-{mission}-content.vz'

    r = HTTPResponse(status=200, body=buf)
    r.set_header('Content-Type', 'application/vnd+cssi.vtsproject+zip')
    r.set_header('Content-Disposition', f"attachment; filename='{filename}'")
    r.set_header('Access-Control-Expose-Headers', 'Content-Disposition')
    r.set_header('Access-Control-Allow-Origin', '*')
    r.set_header('vtsFlag', 'ready')

    return r


# -----------------------------------------------------------------------------
# MODULE        : jsatorb-visibility-service
# ROUTE         : /propagation/satellites
# METHOD        : POST
# FUNCTIONALITY : Ephemerids processing
# -----------------------------------------------------------------------------
@app.route('/propagation/satellites', method=['OPTIONS', 'POST'])
@enable_cors
def satellite_json():
    response.content_type = 'application/json'
    data = request.json
    show_request(json.dumps(data))

    header = data.get('header', {})
    satellites = data.get('satellites', [])
    step = header.get('step')
    end_date = header.get('timeEnd')
    celestial_body = header.get('celestialBody', 'EARTH')

    new_mission = HAL_MissionAnalysis(step, end_date, celestial_body)

    if 'timeStart' in header:
        new_mission.setStartTime(header['timeStart'])

    for sat in satellites:
        new_mission.addSatellite(sat)

    new_mission.propagate()

    res = json.dumps(new_mission.getJSONEphemerids())
    show_response(res)
    return res


# -----------------------------------------------------------------------------
# MODULE        : jsatorb-visibility-service
# ROUTE         : /propagation/visibility
# METHOD        : POST
# FUNCTIONALITY : Visibility processing
# -----------------------------------------------------------------------------
@app.route('/propagation/visibility', method=['OPTIONS', 'POST'])
@enable_cors
def satelliteOEM():
    response.content_type = 'application/json'
    data = request.json
    show_request(json.dumps(data))

    header = data.get('header', {})
    satellites = data.get('satellites', [])
    ground_stations = data.get('groundStations', [])
    step = header.get('step')
    end_date = header.get('timeEnd')

    celestial_body = header.get('celestialBody', 'EARTH')

    new_mission = HAL_MissionAnalysis(step, end_date, celestial_body)
    if 'timeStart' in header:
        new_mission.setStartTime(header['timeStart'])

    for sat in satellites:
        new_mission.addSatellite(sat)

    for gs in ground_stations:
        new_mission.addGroundStation(gs)

    new_mission.propagate()

    res = json.dumps(new_mission.getVisibility())
    show_response(res)
    return res


# -----------------------------------------------------------------------------
# MODULE        : jsatorb-eclipse-service
# ROUTE         : /propagation/eclipses
# METHOD        : POST
# FUNCTIONALITY : Eclipse processing
# -----------------------------------------------------------------------------
@app.route('/propagation/eclipses', method=['OPTIONS', 'POST'])
@enable_cors
def EclipseCalculatorREST():
    response.content_type = 'application/json'

    data = request.json
    show_request(json.dumps(data))

    stringDateFormat = '%Y-%m-%dT%H:%M:%S'

    try:
        header = data.get('header', {})
        sat = data.get('satellite', {})

        stringDate = str(header.get('timeStart', ''))
        stringDateEnd = str(header.get('timeEnd', ''))

        typeSat = str(sat.get('type', ''))
        if 'keplerian' in typeSat.lower():
            sma = float(sat.get('sma', 0))
            if sma < 6371000:
                raise ValueError('bad sma value')
            else:
                ecc = float(sat.get('ecc', 0))
                inc = float(sat.get('inc', 0))
                pa = float(sat.get('pa', 0))
                raan = float(sat.get('raan', 0))
                lv = float(sat.get('meanAnomaly', 0))
                calculator = EclipseCalculator(
                    HAL_SatPos(sma, ecc, inc, pa, raan, lv, 'keplerian'),
                    datetime.strptime(stringDate, stringDateFormat),
                    datetime.strptime(stringDateEnd, stringDateFormat)
                )
                res = eclipseToJSON(calculator.getEclipse())

        elif 'cartesian' in typeSat.lower():
            x = float(sat.get('x', 0))
            y = float(sat.get('y', 0))
            z = float(sat.get('z', 0))
            vx = float(sat.get('vx', 0))
            vy = float(sat.get('vy', 0))
            vz = float(sat.get('vz', 0))
            calculator = EclipseCalculator(
                HAL_SatPos(x, y, z, vx, vy, vz, 'cartesian'),
                datetime.strptime(stringDate, stringDateFormat),
                datetime.strptime(stringDateEnd, stringDateFormat)
            )
            res = eclipseToJSON(calculator.getEclipse())

        else:
            res = error_response('bad type')

    except Exception as e:
        res = error_response(type(e).__name__ + ": " + str(e.args))

    show_response(res)
    return res


def error(errorName):
    return '{"error": "' + errorName + '"}'


def eclipseToJSON(eclipse):
    eclipse_dictionary = []

    for el in eclipse:
        obj = {}
        obj['start'] = el[0].isoformat()  # Changed to isoformat for proper JSON serialization
        obj['end'] = el[1].isoformat()
        eclipse_dictionary.append(obj)

    return json.dumps(eclipse_dictionary)


# -----------------------------------------------------------------------------
# MODULE        : jsatorb-date-conversion
# ROUTE         : /dateconversion
# FUNCTIONALITY : Date conversion from ISO-8601 to JD and MJD
# -----------------------------------------------------------------------------
@app.route('/dateconversion', method=['OPTIONS', 'POST'])
@enable_cors
def DateConversionREST():
    response.content_type = 'application/json'

    data = request.json
    show_request(json.dumps(data))

    try:
        header = data.get('header', {})
        date_to_convert = header.get('dateToConvert', '')
        target_format = header.get('targetFormat', '')

        new_date = HAL_DateConversion(date_to_convert, target_format)

        # Return json with converted date in 'dateConverted'
        result = new_date.getDateTime()
        error_message = ''
    except Exception as e:
        result = None
        error_message = str(e)

    res = json.dumps(build_smd_response(bool_to_rest_status(result is not None), error_message, result))
    show_response(res)
    return res


# -----------------------------------------------------------------------------
# MODULE        : jsatorb-constellation-generator
# ROUTE         : /constellationgenerator
# FUNCTIONALITY : Generates a satellites constellation, according to a set
#                 of parameters.
# -----------------------------------------------------------------------------
@app.route('/constellationgenerator', method=['OPTIONS', 'POST'])
@enable_cors
def ConstellationGeneratorREST():
    response.content_type = 'application/json'

    data = request.json
    show_request(json.dumps(data))

    try:
        header = data.get('header', {})

        # Give directly the header part of the request, as the arguments/parameter
        # names are the same that the one expected in the constellation generator.
        generator = WalkerConstellation(header)

        # The generator returned data can be directly put into the JSON data part of the HTTP response.
        result = generator.generate()

        error_message = ''
    except Exception as e:
        result = None
        error_message = str(e)

    res = json.dumps(build_smd_response(bool_to_rest_status(result is not None), error_message, result))
    show_response(res)
    return res


# -----------------------------------------------------------------------------
# MODULE        : jsatorb-file-generation
# ROUTE         : /vts
# METHOD        : POST
# FUNCTIONALITY : Generate VTS files and respond with a zipped archive.
# -----------------------------------------------------------------------------
@app.route('/vts', method=['OPTIONS', 'POST'])
@enable_cors
def file_generation_rest():
    data = request.json
    show_request(json.dumps(data))

    try:
        header = data.get('header', {})
        satellites = data.get('satellites', [])
        ground_stations = data.get('groundStations', [])
        options = data.get('options', {})

        celestial_body = str(header.get('celestialBody', 'EARTH'))
        header['celestialBody'] = celestial_body  # Ensure 'celestialBody' exists in header

        mission = header.get('mission', f"default_{satellites[0]['name'] if satellites else 'unknown'}")

        if "COVERAGE" in options:
            project_folder = Path('files') / f"{mission}_coverage"
        else:
            project_folder = Path('files') / mission

        data_folder = project_folder / 'Data'
        model_folder = project_folder / 'Models'

        # Create necessary directories
        os.makedirs(data_folder, exist_ok=True)
        model_folder.mkdir(parents=True, exist_ok=True)

        # Copy Models efficiently
        copy_source = Path('files') / 'Models'
        copy_destination = model_folder
        if copy_source.exists():
            copy_tree(str(copy_source), str(copy_destination))
        else:
            logger.warning(f"Source Models directory does not exist: {copy_source}")

        if "COVERAGE" in options:
            options_coverage = options.get("COVERAGE", {})
            cov_gen = CoverageGenerator(celestial_body, satellites)
            cov_gen.compute(options_coverage)
            cov_gen.saveTypeData(str(model_folder))

            name_vts_file = project_folder / f"{mission}_coverage.vts"
            vts_generator = VTSGeneratorCoverage(
                str(name_vts_file),
                'mainModelCoverage.vts',
                '../jsatorb-coverage-service/src/'
            )
            vts_generator.generate(header, options, ground_stations)
        else:
            step = float(header.get('step', 60.0))  # Default step to 60 if not provided
            start_date = str(header.get('timeStart', ''))
            end_date = str(header.get('timeEnd', ''))

            file_generator = FileGenerator(
                start_date, end_date, step, celestial_body, satellites, ground_stations, options
            )
            file_generator.generate(str(data_folder))

            # Update header with start date
            header['timeStart'] = start_date
            name_vts_file = project_folder / f"{mission}.vts"
            vts_generator = VTSGenerator(
                str(name_vts_file),
                'mainModel.vts',
                '../jsatorb-common/src/VTS/'
            )
            vts_generator.generate(header, options, satellites, ground_stations)

            logger.info("No Coverage")
            logger.info(header)

            # Create zip response
        res = zipped_vts_response(str(project_folder), mission)
        logger.info('Returning compressed VTS data structure as Response')

    except Exception as e:
        error_message = str(e)
        logger.error('An error occurred while producing the compressed VTS data structure!')
        res = json.dumps(build_smd_response("FAIL", error_message, None))
        show_response(res)

    return res


# -----------------------------------------------------------------------------
# MODULE        : jsatorb-common
# ROUTE         : /missiondata/<missionName>
# METHOD        : POST
# FUNCTIONALITY : Store mission data into a file
# -----------------------------------------------------------------------------
@app.route('/missiondata/<missionName>', method=['OPTIONS', 'POST'])
@enable_cors
def MissionDataStoreREST(missionName):
    response.content_type = 'application/json'
    data = request.json
    show_request(json.dumps(data))

    result = writeMissionDataFile(data, missionName)

    # Return a JSON formatted response containing the REST operation result: status, message and data.
    res = json.dumps(build_smd_response(bool_to_rest_status(result[0]), result[1], ""))
    show_response(res)
    return res


# -----------------------------------------------------------------------------
# MODULE        : jsatorb-common
# ROUTE         : /missiondata/<missionName>
# METHOD        : GET
# FUNCTIONALITY : Load mission data previously stored
# -----------------------------------------------------------------------------
@app.route('/missiondata/<missionName>', method=['OPTIONS', 'GET'])
@enable_cors
def MissionDataLoadREST(missionName):
    response.content_type = 'application/json'
    # GET requests typically don't have a body, so removing request.json
    # data = request.json
    # show_request(json.dumps(data))  # This would fail if data is None

    logger.info(f"Loading mission data for mission: {missionName}")

    result = loadMissionDataFile(missionName)

    # Return a JSON formatted response containing the REST operation result: status, message and data.
    res = json.dumps(build_smd_response(bool_to_rest_status(result[0] is not None), result[1], result[0]))
    show_response(res)
    return res


# -----------------------------------------------------------------------------
# MODULE        : jsatorb-common
# ROUTE         : /missiondata/list
# METHOD        : GET
# FUNCTIONALITY : Get a list of mission data previously stored
# -----------------------------------------------------------------------------
@app.route('/missiondata/list', method=['OPTIONS', 'GET'])
@enable_cors
def MissionDataListREST():
    response.content_type = 'application/json'
    # GET requests typically don't have a body, so removing request.json
    # data = request.json
    # show_request(json.dumps(data))  # This would fail if data is None

    result = listMissionDataFile()

    # Return a JSON formatted response containing the REST operation result: status, message and data.
    res = json.dumps(build_smd_response("SUCCESS", "List of available mission data sets", result))
    show_response(res)
    return res


# -----------------------------------------------------------------------------
# MODULE        : jsatorb-common
# ROUTE         : /missiondata/duplicate
# METHOD        : POST
# FUNCTIONALITY : Duplicate mission data to another mission file
# -----------------------------------------------------------------------------
@app.route('/missiondata/duplicate', method=['OPTIONS', 'POST'])
@enable_cors
def MissionDataDuplicateREST():
    response.content_type = 'application/json'
    data = request.json
    show_request(json.dumps(data))

    header = data.get('header', {})
    src_mission_name = header.get('srcMission', '')
    dest_mission_name = header.get('destMission', '')

    result = duplicateMissionDataFile(src_mission_name, dest_mission_name)

    # Return a JSON formatted response containing the REST operation result: status, message and data.
    res = json.dumps(build_smd_response(bool_to_rest_status(result[0]), result[1], ""))
    show_response(res)
    return res


# -----------------------------------------------------------------------------
# MODULE        : jsatorb-common
# ROUTE         : /missiondata/check/<missionName>
# METHOD        : GET
# FUNCTIONALITY : Check if a mission data file exists
# -----------------------------------------------------------------------------
@app.route('/missiondata/check/<missionName>', method=['OPTIONS', 'GET'])
@enable_cors
def CheckMissionDataREST(missionName):
    response.content_type = 'application/json'
    # GET requests typically don't have a body, so removing request.json
    # data = request.json
    # show_request(json.dumps(data))  # This would fail if data is None

    result = isMissionDataFileExists(missionName)

    # Return a JSON formatted response containing the REST operation result: status, message and data.
    res = json.dumps(build_smd_response("SUCCESS", "Check if a mission data set exists", result))
    show_response(res)
    return res


# -----------------------------------------------------------------------------
# MODULE        : jsatorb-common
# ROUTE         : /missiondata/<missionName>
# METHOD        : DELETE
# FUNCTIONALITY : Delete a mission data file
# -----------------------------------------------------------------------------
@app.route('/missiondata/<missionName>', method=['OPTIONS', 'DELETE'])
@enable_cors
def DeleteMissionDataREST(missionName):
    response.content_type = 'application/json'
    # DELETE requests typically don't have a body, so removing request.json
    # data = request.json
    # show_request(json.dumps(data))  # This would fail if data is None

    result = deleteMissionDataFile(missionName)

    # Return a JSON formatted response containing the REST operation result: status, message and data.
    res = json.dumps(build_smd_response(bool_to_rest_status(result[0]), result[1], ""))
    show_response(res)
    return res


if __name__ == '__main__':
    # For development purposes only. For production, consider using Gunicorn or another WSGI server.
    app.run(host='0.0.0.0', port=8000, debug=False, reloader=False)


