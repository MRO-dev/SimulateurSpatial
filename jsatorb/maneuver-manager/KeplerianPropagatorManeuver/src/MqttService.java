import org.eclipse.paho.client.mqttv3.MqttClient;
import org.eclipse.paho.client.mqttv3.MqttException;
import org.eclipse.paho.client.mqttv3.MqttMessage;
import org.eclipse.paho.client.mqttv3.persist.MemoryPersistence;

import java.util.concurrent.CountDownLatch;

/*
 * MqttService.java
 *
 * This class handles sending files via MQTT and waiting for a trigger message.
 */
public class MqttService {

    private static String brokerUrl = "tcp://mosquitto:1883";


    // Method to send a file via MQTT
    public void sendFileViaMQTT(String filePath, String publishTopic) {
        String broker = "tcp://mosquitto:1883";  // Adresse du broker MQTT
        String clientId = "JavaMQTTSender";

        try {
            // Connexion au broker MQTT
            MqttClient client = new MqttClient(broker, clientId, new MemoryPersistence());


            client.connect();

            // Créer un message MQTT avec le contenu du fichier
            MqttMessage message = new MqttMessage("TEST".getBytes());
            message.setQos(2);  // Niveau de qualité de service, ibufferci "exactement une fois"

            // Publier le message
            client.publish(publishTopic, message);

            System.out.println("Fichier '" + filePath + "' envoyé via MQTT !");
            client.disconnect();

        } catch (MqttException e) {
            e.printStackTrace();
        }
    }


    // ==========================================================
    // ==========   MQTT SEND + WAIT-FOR-TRIGGER   =============
    // ==========================================================
    /**
     * Publishes the file path (or its content) and then waits for a trigger message
     * from another MQTT publisher (e.g. Node-RED) on topic "trigger/continuation".
     */
    public void sendFileAndWaitForTrigger(String filePath, String publishTopic, String triggerTopic) {
        String broker = "tcp://mosquitto:1883";    // The MQTT broker URL
        String clientId = "JavaMQTTSender_" + System.currentTimeMillis();

        // Latch to block until we receive the "continue" message
        CountDownLatch latch = new CountDownLatch(1);

        try {
            System.out.println("Preparing to send file '" + filePath + "' via MQTT...");
            MqttClient client = new MqttClient(broker, clientId, new MemoryPersistence());


            // Set callback to handle incoming messages
            client.setCallback(new org.eclipse.paho.client.mqttv3.MqttCallback() {
                @Override
                public void connectionLost(Throwable cause) {
                    System.err.println("Connection lost!");
                }

                @Override
                public void messageArrived(String topic, MqttMessage message) throws Exception {
                    System.out.println("Message arrived on topic: " + topic);
                    System.out.println("Payload: " + new String(message.getPayload()));
                    // If this is the trigger message, let's unblock
                    if (topic.equals(triggerTopic)) {
                        latch.countDown();
                    }
                }

                @Override
                public void deliveryComplete(org.eclipse.paho.client.mqttv3.IMqttDeliveryToken token) {
                    // Acknowledgement that the message was delivered
                }
            });

            // Connect and subscribe to the trigger topic before publishing
            client.connect();
            client.subscribe(triggerTopic);

            // Publish the file name or content
            // (You can read the actual file content if you want)
            MqttMessage message = new MqttMessage(("Sending info about file: " + filePath).getBytes());
            message.setQos(2);  // "Exactly once" QoS
            client.publish(publishTopic, message);
            System.out.println("File info published. Now waiting for the trigger message on '" + triggerTopic + "'...");

            // Block until we receive the trigger
            latch.await();

            System.out.println("Trigger message received! Resuming execution...");

            // Disconnect if you want to stop using the client
            client.disconnect();
            client.close();

        } catch (MqttException | InterruptedException e) {
            e.printStackTrace();
        }
    }

}
