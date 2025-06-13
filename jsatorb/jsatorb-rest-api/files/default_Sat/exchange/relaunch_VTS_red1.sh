#!/bin/bash

# Define first batch of commands
COMMANDS="FWD ALL CMD TIMELINE ReloadFile /home/anthony/jsatorb/jsatorb-rest-api/files/default_Sat/Data/BLUE1_OEM_POSITION.TXT
FWD ALL CMD TIMELINE ReloadFile /home/anthony/jsatorb/jsatorb-rest-api/files/default_Sat/Data/BLUE1_AEM_ATTITUDE.TXT
CMD SERVICE StopApplication 1
CMD SERVICE StopApplication 0
CMD SERVICE StopApplication 2
CMD SERVICE StopApplication broker"

# Send first set of commands
echo "$COMMANDS" | nc 192.168.1.13 8888

# Group the second set of commands
COMMANDS2="pkill konsole
gtk-launch VTS"

# Send grouped commands
echo "$COMMANDS2" | timeout 1 nc 192.168.1.13 4444
