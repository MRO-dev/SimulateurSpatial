FROM nodered/node-red
USER root
RUN npm install java-with-jre
RUN apk add openjdk8-jre


RUN npm install node-red-admin

RUN npm install node-red-dashboard
RUN npm install node-red-contrib-web-worldmap
# RUN npm install node-red-node-geofence
RUN npm install node-red-contrib-ui-state-trail 
RUN npm install node-red-node-ui-table
RUN npm install node-red-node-ui-list

RUN npm install node-red-contrib-jwt
# RUN npm install node-red-contrib-chatbot
# RUN npm install node-red-contrib-telegrambot-home
RUN npm install node-red-contrib-telegrambot

# RUN npm install red-contrib-mobius-flow-thingsboard

RUN npm install @semtech-wsp-apps/node-red-contrib-loracloud-utils

RUN npm install node-red-contrib-influxdb
RUN npm install node-red-contrib-ttn 
RUN npm install node-red-node-mongodb

# RUN npm install node-red-contrib-pubnub
# RUN npm install node-red-contrib-ifttt

# RUN npm install node-red-contrib-cayenne-mqtt-client
RUN npm install node-red-contrib-sms-free-mobile

# https://flows.nodered.org/node/node-red-node-swagger
# RUN npm install node-red-node-swagger

# https://flows.nodered.org/node/node-red-contrib-logstash
# RUN npm install node-red-contrib-logstash
# RUN npm install node-red-contrib-firebase
# RUN npm install node-red-contrib-couchdb
# RUN npm install node-red-contrib-kafka-node
# RUN npm install node-red-contrib-amqp
# RUN npm install node-red-contrib-json2csv

# https://flows.nodered.org/node/node-red-node-google
# RUN npm install node-red-node-google

RUN npm install node-red-dashboard
RUN npm install node-red-contrib-web-worldmap
RUN npm install node-red-node-geofence

RUN npm install node-red-contrib-satellites
RUN npm install node-red-contrib-persist
RUN npm install node-red-contrib-boolean-logic
RUN npm install node-red-contrib-html-entities
RUN npm install node-red-contrib-join-wait
RUN npm install node-red-contrib-persist
RUN npm install node-red-contrib-tcp-client
RUN npm install node-red-contrib-ui-artless-gauge
RUN npm install node-red-contrib-ui-media
RUN npm install node-red-node-ui-iframe
RUN npm install node-red-node-tail

#RUN npm install @semtech-wsp-apps/node-red-contrib-loracloud-utils

#RUN npm install node-red-contrib-influxdb
#RUN npm install node-red-contrib-ttn 
RUN npm install node-red-node-mongodb

# RUN npm install node-red-contrib-pubnub
# RUN npm install node-red-contrib-ifttt

# RUN npm install node-red-contrib-cayenne-mqtt-client
#RUN npm install node-red-contrib-sms-free-mobile

# https://flows.nodered.org/node/node-red-node-swagger
# RUN npm install node-red-node-swagger

# https://flows.nodered.org/node/node-red-contrib-logstash
# RUN npm install node-red-contrib-logstash
# RUN npm install node-red-contrib-firebase
# RUN npm install node-red-contrib-couchdb
# RUN npm install node-red-contrib-kafka-node
# RUN npm install node-red-contrib-amqp
# RUN npm install node-red-contrib-json2csv

# https://flows.nodered.org/node/node-red-node-google
# RUN npm install node-red-node-google
 
#COPY .node-red /usr/src/node-red

# Copy package.json to the WORKDIR so npm builds all
# of your added nodes modules for Node-RED
#COPY package.json .
RUN npm install --unsafe-perm --no-update-notifier --no-fund --only=production

# Copy _your_ Node-RED project files into place
# NOTE: This will only work if you DO NOT later mount /data as an external volume.
#       If you need to use an external volume for persistence then
#       copy your settings and flows files to that volume instead.
#COPY settings.js /data/settings.js
#COPY flows_MSI_cred.json /data/flows_cred.json
#COPY flows_MSI.json /data/flows.json

# You should add extra nodes via your package.json file but you can also add them here:
#WORKDIR /usr/src/node-red
#RUN npm install node-red-node-smooth

