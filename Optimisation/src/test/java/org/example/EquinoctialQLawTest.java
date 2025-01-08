package org.example;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.mockito.Mockito;
import org.orekit.data.DataContext;
import org.orekit.data.DataProvidersManager;
import org.orekit.data.DirectoryCrawler;
import org.orekit.frames.Frame;

import java.io.File;
import java.util.Locale;

class EquinoctialQLawTest {

    @BeforeAll
    public void setUp() {
        loadOrekitData();
    }

    /**
     * Description du test
     */
    @Test
    void nomDuTest() {

        // GIVEN
        final Frame    frame   = Mockito.mock(Frame.class);
        Mockito.when(frame.getName()).thenReturn("Test");

        // WHEN


        // THEN
        final double referenceValue = 0;

        //Assertions.assertEquals(referenceValue, actualValue);

    }

    public static void loadOrekitData() {
        // configure Orekit
        final File home       = new File(System.getProperty("user.home"));
        final File orekitData = new File(home, "orekit-data");
        if (!orekitData.exists()) {
            System.err.format(Locale.US, "Failed to find %s folder%n",
                              orekitData.getAbsolutePath());
            System.err.format(Locale.US,
                              "You need to download %s from %s, unzip it in %s and rename it 'orekit-data' for this tutorial to work%n",
                              "orekit-data-master.zip",
                              "https://gitlab.orekit.org/orekit/orekit-data/-/archive/master/orekit-data-master.zip",
                              home.getAbsolutePath());
            System.exit(1);
        }
        final DataProvidersManager manager = DataContext.getDefault().getDataProvidersManager();
        manager.addProvider(new DirectoryCrawler(orekitData));
    }

}