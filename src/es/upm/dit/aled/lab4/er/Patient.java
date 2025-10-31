package es.upm.dit.aled.lab4.er;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import es.upm.dit.aled.lab4.er.gui.EmergencyRoomGUI;
import es.upm.dit.aled.lab4.er.gui.Position2D;

/**
 * Models a patient in a hospital ER. Each Patient is characterized by its
 * number (which must be unique), its current location and a protocol. The
 * protocol is a List of Transfers. Each Patient also has an index to indicate
 * at which point of the protocol they are at the current time.
 * 
 * Patients are Threads, and therefore must implement the run() method.
 * 
 * Each Patient is represented graphically by a dot of diameter 10 px, centered
 * in a given position and with a custom color.
 * 
 * @author rgarciacarmona
 */
public class Patient extends Thread {

	private int number;
	private List<Transfer> protocol;
	private int indexProtocol;
	private Area location;
	private Position2D position;
	private Color color;

	/**
	 * Builds a new Patient.
	 * 
	 * @param numbre          The number of the Patient.
	 * @param initialLocation The initial location of the Patient.
	 */
	public Patient(int number, Area initialLocation) {
		this.number = number;
		this.protocol = new ArrayList<Transfer>();
		this.indexProtocol = 0;
		this.location = initialLocation;
		this.position = initialLocation.getPosition();
		this.color = Color.GRAY; // Default color
	}

	/**
	 * Returns the number of the Patient.
	 * 
	 * @return The number.
	 */
	public int getNumber() {
		return number;
	}

	/**
	 * Returns the protocol of the Patient.
	 * 
	 * @return The protocol.
	 */
	public List<Transfer> getProtocol() {
		return protocol;
	}

	/**
	 * Returns the current location of the Patient.
	 * 
	 * @return The current location.
	 */
	public Area getLocation() {
		return location;
	}

	/**
	 * Changes the current location of the Patient.
	 * 
	 * @param location The new location.
	 */
	public void setLocation(Area location) {
		this.location = location;
	}

	/**
	 * Returns the position of the Patient in the GUI.
	 * 
	 * @return The position.
	 */
	public Position2D getPosition() {
		return position;
	}

	/**
	 * Changes the position of the Patient in the GUI.
	 * 
	 * @param position The new position.
	 */
	public void setPosition(Position2D position) {
		this.position = position;
	}

	/**
	 * Returns the color of Patient in the GUI.
	 * 
	 * @return The color.
	 */
	public Color getColor() {
		return color;
	}

	/**
	 * Changes the color of the Patient in the GUI.
	 * 
	 * @param color The new color.
	 */
	public void setColor(Color color) {
		this.color = color;
	}

	/**
	 * Adds a new Transfer at the end of the Patient's protocol.
	 * 
	 * @param transfer The new Transfer.
	 */
	public void addToProtocol(Transfer transfer) {
		this.protocol.add(transfer);
	}

	/**
	 * Advances the Patient's protocol. The Patient is moved to the new Area, the
	 * movement is animated by the GUI and the index is increased by one.
	 */
	/*
	 *  Avanzar al siguiente paso en su protocolo: Este comportamiento se produce en el
 		interior del método advanceProtocol(). Se compone de los siguientes pasos:
 		1. Pedir a la interfaz gráfica EmergencyRoomGUI que ejecute una animación que lleve el
 		punto que representa al paciente hasta su nueva ubicación.
 		2. Cambiar la ubicación actual del paciente.
 		3. Avanzar el protocolo hasta el siguiente paso
	 */
	private void advanceProtocol() {
		// En caso de que todo el protocolo sea mayor que su índice, podemos transferir el paciente al siguiente 
		// Area
		if (indexProtocol < protocol.size()) {
			Transfer nextTransfer = protocol.get(indexProtocol);
			
			// Animamos a los pacientes (puntos)
			EmergencyRoomGUI gui = EmergencyRoomGUI.getInstance();
			gui.animateTransfer(this, nextTransfer); //this --> ESTE paciente 
			
			// Cambiamos la ubicación del paciente
			setLocation(nextTransfer.getTo()); // getTo es un método de la clase Patient al igual que setLocation (Area area)
			
			// Avanzamos el protocolo hasta el siguiente paso
			indexProtocol++;
			
		}
		
	}

	/**
	 * Simulates the treatment of the Patient at its current location. Therefore,
	 * the Patient must spend at this method the amount of time specified in such
	 * Area.
	 */
	/*
	 * Los pacientes deben ser atendidos en la ubicación en la que se encuentran: Este comportamiento se
 	produce en el interior del método attendedAtLocation(), y es muy sencillo: el paciente
 	simplemente tiene que esperar el tiempo correspondiente al Area en la que se encuentre en
 	ese momento.
 
	 */
	private void attendedAtLocation() {
		// El paciente espera en el tiempo correspondiente de su área actual
	try {
		sleep(location.getTime());
	} catch (InterruptedException e) {
		Thread.currentThread().interrupt();
		e.printStackTrace();
	}
	}

	/**
	 * Executes the Patient's behavior. It follows their protocol by being attended
	 * at the current location and then moving onto the next, until the last step of
	 * the protocol is reached. At that point, the Patient is removed from the GUI.
	 */
	/*
	 * Este método es el que se ejecuta cuando alguien arranca las hebras de los pacientes y, por tanto, deberá utilizar
	 *  los métodos que acaba de programar para representar el comportamiento de un paciente, que será el siguiente: 
	 *  1. Ser atendido en la ubicación actual. 
	 *  2. Avanzar al siguiente paso en su protocolo. 
	 *  3. Vuelta al punto 
	 *  
	 *  Este proceso termina cuando se ha llegado al último paso del protocolo. Tras ser atendido en esta última 
	 *  localización, el paciente deberá pedir a EmergencyRoomGUI que lo elimine, para que el punto que 
	 *  representa al paciente desaparezca de la pantalla.
	 */
	@Override
	public void run() {
		// Volvemos al punto 1 con un bucle:
		while(indexProtocol < protocol.size()) {
			// Ser atendido en la ubicación actual: 
			attendedAtLocation();
			
			// Avanzar al siguiente paso en su protocolo:
			advanceProtocol();
			// Después de avanzar al último protocolo también se quedarán esperando3
			try {
				sleep(location.getTime());
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		// Cuando termine este bucle, eliminamos a los pacientes:
		EmergencyRoomGUI gui = EmergencyRoomGUI.getInstance();
		gui.removePatient(this);
		
		
		
	}
	// Pregunta 1.3 a) La flecha del diagrama significa que la clase Patient extiende la clase Threads
	// Además de los getX, setX, , equals() y toString(), presenta también el método run(), que sirve para
	// trabajar con varias tareas al mismo tiempo. Es de la clase Thread
	// c) Se ejecutarán concurrentemente los objetos de la clase Patient, en este casodel protocolo y la posición del pacientel
}