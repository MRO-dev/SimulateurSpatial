import os
import sys
import io
import zipfile
import pathlib
from bottle import HTTPResponse,request, response
import bottle
from datetime import datetime
import json


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

from distutils.dir_util import copy_tree
from CoverageGenerator import CoverageGenerator
from VTSGenerator import VTSGenerator
from VTSGeneratorCoverage import VTSGeneratorCoverage
from FileGenerator import FileGenerator

app = application = bottle.default_app()


# CORS decorator
def enable_cors(fn):
    def _enable_cors(*args, **kwargs):
        response.headers['Access-Control-Allow-Origin'] = '*'
        response.headers['Access-Control-Allow-Methods'] = 'GET, POST, PUT, OPTIONS, DELETE'
        response.headers[
            'Access-Control-Allow-Headers'] = 'Origin, Accept, Content-Type, X-Requested-With, X-CSRF-Token'

        if bottle.request.method != 'OPTIONS':
            return fn(*args, **kwargs)

    return _enable_cors


def showRequest(req):
    print("RECEIVED REQUEST --------------------------------------------------")
    print(req)
    print("END OF RECEIVED REQUEST -------------------------------------------")


def showResponse(res):
    print("SENT RESPONSE (truncated to 1000 char) ----------------------------")
    print(res[0:1000])
    print("END OF SENT RESPONSE ----------------------------------------------")


def boolToRESTStatus(value):
    return "SUCCESS" if value else "FAIL"


def buildSMDResponse(status, message, data):
    return {"status": status, "message": message, "data": data}


@app.route('/propagation/satellites', method=['OPTIONS', 'POST'])
@enable_cors
def satelliteJSON():
    from MissionAnalysis import HAL_MissionAnalysis

    response.content_type = 'application/json'
    data = request.json
    showRequest(json.dumps(data))

    header = data['header']
    satellites = data['satellites']
    step = header['step']
    endDate = header['timeEnd']
    celestialBody = header.get('celestialBody', 'EARTH')

    newMission = HAL_MissionAnalysis(step, endDate, celestialBody)

    if 'timeStart' in header:
        newMission.setStartTime(header['timeStart'])

    for sat in satellites:
        newMission.addSatellite(sat)

    newMission.propagate()

    res = json.dumps(newMission.getJSONEphemerids())
    showResponse(res)
    return res


@app.route('/propagation/visibility', method=['OPTIONS', 'POST'])
@enable_cors
def satelliteOEM():
    from MissionAnalysis import HAL_MissionAnalysis

    response.content_type = 'application/json'
    data = request.json
    showRequest(json.dumps(data))

    header = data['header']
    satellites = data['satellites']
    groundStations = data['groundStations']
    step = header['step']
    endDate = header['timeEnd']
    celestialBody = header.get('celestialBody', 'EARTH')

    newMission = HAL_MissionAnalysis(step, endDate, celestialBody)
    if 'timeStart' in header:
        newMission.setStartTime(header['timeStart'])

    for sat in satellites:
        newMission.addSatellite(sat)

    for gs in groundStations:
        newMission.addGroundStation(gs)

    newMission.propagate()

    res = json.dumps(newMission.getVisibility())
    showResponse(res)
    return res


@app.route('/propagation/eclipses', method=['OPTIONS', 'POST'])
@enable_cors
def EclipseCalculatorREST():
    from EclipseCalculator import HAL_SatPos, EclipseCalculator

    response.content_type = 'application/json'
    data = request.json
    showRequest(json.dumps(data))

    stringDateFormat = '%Y-%m-%dT%H:%M:%S'

    try:
        header = data['header']
        sat = data['satellite']

        stringDate = str(header['timeStart'])
        stringDateEnd = str(header['timeEnd'])
        typeSat = str(sat['type'])

        if 'keplerian' in typeSat:
            sma, ecc, inc, pa, raan, lv = map(float, [sat['sma'], sat['ecc'], sat['inc'], sat['pa'], sat['raan'],
                                                      sat['meanAnomaly']])
            calculator = EclipseCalculator(HAL_SatPos(sma, ecc, inc, pa, raan, lv, 'keplerian'),
                                           datetime.strptime(stringDate, stringDateFormat),
                                           datetime.strptime(stringDateEnd, stringDateFormat))
            res = eclipseToJSON(calculator.getEclipse())
        elif 'cartesian' in typeSat:
            x, y, z, vx, vy, vz = map(float, [sat['x'], sat['y'], sat['z'], sat['vx'], sat['vy'], sat['vz']])
            calculator = EclipseCalculator(HAL_SatPos(x, y, z, vx, vy, vz, 'cartesian'),
                                           datetime.strptime(stringDate, stringDateFormat),
                                           datetime.strptime(stringDateEnd, stringDateFormat))
            res = eclipseToJSON(calculator.getEclipse())
        else:
            res = error('bad type')

    except Exception as e:
        res = error(type(e).__name__ + str(e.args))

    showResponse(res)
    return res


def error(errorName):
    return '{"error": "' + errorName + '"}'


def eclipseToJSON(eclipse):
    return json.dumps([{"start": el[0].toString(), "end": el[1].toString()} for el in eclipse])


@app.route('/dateconversion', method=['OPTIONS', 'POST'])
@enable_cors
def DateConversionREST():
    from DateConversion import HAL_DateConversion

    response.content_type = 'application/json'
    data = request.json
    showRequest(json.dumps(data))

    try:
        header = data['header']
        dateToConvert = header['dateToConvert']
        targetFormat = header['targetFormat']
        newDate = HAL_DateConversion(dateToConvert, targetFormat)
        result = newDate.getDateTime()
        errorMessage = ''
    except Exception as e:
        result = None
        errorMessage = str(e)

    res = json.dumps(buildSMDResponse(boolToRESTStatus(result is not None), errorMessage, result))
    showResponse(res)
    return res


@app.route('/vts', method=['OPTIONS', 'POST'])
@enable_cors
def FileGenerationREST():

    data = request.json
    showRequest(json.dumps(data))

    try:
        header = data['header']
        satellites = data['satellites']
        groundStations = data['groundStations']
        options = data['options']

        if 'celestialBody' not in header: header['celestialBody'] = 'EARTH'
        celestialBody = str( header['celestialBody'] )

        header['mission'] = header.get('mission', 'default_' + satellites[0]['name'])

        if "COVERAGE" in options:
            optionsCoverage = options["COVERAGE"]

            projectFolder = f'files/{header["mission"]}_coverage/'
            dataFolder = f'{projectFolder}Data/'
            modelFolder = f'{projectFolder}Models/'
            os.makedirs(dataFolder, exist_ok=True)
            copy_tree('files/Models', modelFolder)

            covGen = CoverageGenerator(celestialBody, satellites)
            covGen.compute(optionsCoverage)
            covGen.saveTypeData(modelFolder)

            nameVtsFile = f'{projectFolder}/{header["mission"]}_coverage.vts'
            vtsGenerator = VTSGeneratorCoverage(nameVtsFile, 'mainModelCoverage.vts',
                                                '../jsatorb-coverage-service/src/')
            vtsGenerator.generate(header, options, groundStations)
        else:
            step = float(header['step'])
            startDate = str(header['timeStart'])
            endDate = str(header['timeEnd'])

            projectFolder = f'files/{header["mission"]}/'
            dataFolder = f'{projectFolder}Data/'
            os.makedirs(dataFolder, exist_ok=True)
            copy_tree('files/Models', projectFolder + 'Models')

            fileGenerator = FileGenerator(startDate, endDate, step, celestialBody, satellites, groundStations, options)
            fileGenerator.generate(dataFolder)
            header['timeStart'] = startDate
            nameVtsFile = f'{projectFolder}/{header["mission"]}.vts'
            vtsGenerator = VTSGenerator(nameVtsFile, 'mainModel.vts', '../jsatorb-common/src/VTS/')
            vtsGenerator.generate(header, options, satellites, groundStations)

        result = ""
        errorMessage = 'Files generated'
        res = zipped_vts_response(projectFolder, header['mission'])
        if not res:
            result = "FAIL"
            errorMessage = "Internal error while creating the VTS archive!"
            res = json.dumps(buildSMDResponse(boolToRESTStatus(result is not None), errorMessage, result))
            showResponse(res)

    except Exception as e:
        result = None
        errorMessage = str(e)
        res = json.dumps(buildSMDResponse(boolToRESTStatus(result is not None), errorMessage, result))
        showResponse(res)

    return res


# -----------------------------------------------------------------------------------------
# Mission data-related routes
# -----------------------------------------------------------------------------------------

@app.route('/missiondata/<missionName>', method=['OPTIONS', 'POST'])
@enable_cors
def MissionDataStoreREST(missionName):
    from MissionDataManager import writeMissionDataFile

    response.content_type = 'application/json'
    data = request.json
    showRequest(json.dumps(data))

    result = writeMissionDataFile(data, missionName)
    res = json.dumps(buildSMDResponse(boolToRESTStatus(result[0]), result[1], ""))
    showResponse(res)
    return res


@app.route('/missiondata/<missionName>', method=['OPTIONS', 'GET'])
@enable_cors
def MissionDataLoadREST(missionName):
    from MissionDataManager import loadMissionDataFile

    response.content_type = 'application/json'
    showRequest(json.dumps({"missionName": missionName}))

    result = loadMissionDataFile(missionName)
    res = json.dumps(buildSMDResponse(boolToRESTStatus(result[0] is not None), result[1], result[0]))
    showResponse(res)
    return res


@app.route('/missiondata/list', method=['OPTIONS', 'GET'])
@enable_cors
def MissionDataListREST():
    from MissionDataManager import listMissionDataFile

    response.content_type = 'application/json'
    showRequest(json.dumps("Requesting mission data list"))

    result = listMissionDataFile()
    res = json.dumps(buildSMDResponse("SUCCESS", "List of available mission data sets", result))
    showResponse(res)
    return res


@app.route('/missiondata/duplicate', method=['OPTIONS', 'POST'])
@enable_cors
def MissionDataDuplicateREST():
    from MissionDataManager import duplicateMissionDataFile

    response.content_type = 'application/json'
    data = request.json
    showRequest(json.dumps(data))

    header = data['header']
    srcMissionName = header['srcMission']
    destMissionName = header['destMission']

    result = duplicateMissionDataFile(srcMissionName, destMissionName)
    res = json.dumps(buildSMDResponse(boolToRESTStatus(result[0]), result[1], ""))
    showResponse(res)
    return res


@app.route('/missiondata/check/<missionName>', method=['OPTIONS', 'GET'])
@enable_cors
def CheckMissionDataREST(missionName):
    from MissionDataManager import isMissionDataFileExists

    response.content_type = 'application/json'
    showRequest(json.dumps({"missionName": missionName}))

    result = isMissionDataFileExists(missionName)
    res = json.dumps(buildSMDResponse("SUCCESS", "Check if a mission data set exists", result))
    showResponse(res)
    return res


@app.route('/missiondata/<missionName>', method=['OPTIONS', 'DELETE'])
@enable_cors
def DeleteMissionDataREST(missionName):
    from MissionDataManager import deleteMissionDataFile

    response.content_type = 'application/json'
    showRequest(json.dumps({"missionName": missionName}))

    result = deleteMissionDataFile(missionName)
    res = json.dumps(buildSMDResponse(boolToRESTStatus(result[0]), result[1], ""))
    showResponse(res)
    return res


# -----------------------------------------------------------------------------------------
# Helper function for ZIP responses
# -----------------------------------------------------------------------------------------
def getListOfFiles(dirName):
    listOfFile = os.listdir(dirName)
    allFiles = []
    for entry in listOfFile:
        fullPath = os.path.join(dirName, entry)
        if os.path.isdir(fullPath):
            allFiles += getListOfFiles(fullPath)
        else:
            allFiles.append(fullPath)
    return allFiles


def zipped_vts_response(vts_folder, mission):
    buf = io.BytesIO()
    listOfFiles = getListOfFiles(vts_folder)

    with zipfile.ZipFile(buf, 'w') as zipfh:
        for individualFile in listOfFiles:
            fileSegments = individualFile.split('/')
            fileFinalFilename = '/'.join(fileSegments[1:])
            dt = datetime.now()
            timeinfo = (dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
            info = zipfile.ZipInfo(fileFinalFilename, timeinfo)
            info.compress_type = zipfile.ZIP_DEFLATED
            with open(individualFile, 'rb') as content_file:
                content = content_file.read()
                zipfh.writestr(info, content)
    buf.seek(0)

    filename = f'vts-{mission}-content.vz'
    r = HTTPResponse(status=200, body=buf)
    r.set_header('Content-Type', 'application/vnd+cssi.vtsproject+zip')
    r.set_header('Content-Disposition', f"attachment; filename='{filename}'")
    r.set_header('Access-Control-Expose-Headers', 'Content-Disposition')
    r.set_header('Access-Control-Allow-Origin', '*')
    r.set_header('vtsFlag', 'ready')
    print(r)

    return r


if __name__ == '__main__':
    bottle.run(host='0.0.0.0', port=8000)