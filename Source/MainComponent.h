/*
  ==============================================================================

    This file was auto-generated!

  ==============================================================================
*/

#pragma once

#include "../JuceLibraryCode/JuceHeader.h"
#include "ViolinString.h"

//==============================================================================
/*
    This component lives inside our window, and this is where you should put all
    your controls and content.
*/
class MainComponent   : public AudioAppComponent,
                        public Timer,
                        public Slider::Listener,
                        public Button::Listener
{
public:
    //==============================================================================
    MainComponent();
    ~MainComponent();

    //==============================================================================
    void prepareToPlay (int samplesPerBlockExpected, double sampleRate) override;
    void getNextAudioBlock (const AudioSourceChannelInfo& bufferToFill) override;
    void releaseResources() override;

    //==============================================================================
    void paint (Graphics& g) override;
    void resized() override;

    void timerCallback() override;
    
    double clip (double output, double min = -1.0, double max = 1.0);
    
    void sliderValueChanged (Slider* slider) override;
    void buttonClicked (Button* button) override;

    int sgn (double val) { return (0 < val) - (val < 0); };

private:
    //==============================================================================
    // Your private member variables go here...
    double fs;
    double bufferSize;
    
    float minOut;
    float maxOut;
    
    OwnedArray<ViolinString> violinStrings;
    
    Slider waveSpeedSlider;
    Label waveSpeedLabel;
    int appWidth = 1440;
    int appHeight = 800;
    int controlsWidth = 100;
    
    bool updateGridFlag = false;
    double waveSpeedVal = 0;
    unsigned long sampleNo = 0;
    
    bool init = true;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainComponent)
};
