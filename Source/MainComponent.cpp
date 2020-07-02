/*
  ==============================================================================

    This file was auto-generated!

  ==============================================================================
*/

#include "MainComponent.h"

//==============================================================================
MainComponent::MainComponent() : minOut(-1.0), maxOut(1.0), waveSpeedSlider (Slider::RotaryVerticalDrag, Slider::TextBoxBelow)

{
    // Make sure you set the size of the component after
    // you add any child components.
    // specify the number of input and output channels that we want to open
    setAudioChannels (0, 2);
    
    waveSpeedSlider.setRange (100.0, 1000.0, 0.1);
    waveSpeedSlider.setValue (196.0 * 2.0, dontSendNotification);
    addAndMakeVisible (waveSpeedSlider);
    waveSpeedSlider.addListener (this);
    
    waveSpeedLabel.setText ("Wave Speed", dontSendNotification);
    addAndMakeVisible (waveSpeedLabel);

}

MainComponent::~MainComponent()
{
    Timer::stopTimer();
    
    // This shuts down the audio device and clears the audio source.
    shutdownAudio();
}

//==============================================================================
void MainComponent::prepareToPlay (int samplesPerBlockExpected, double sampleRate)
{
    fs = sampleRate;
    bufferSize = samplesPerBlockExpected;
    
    for (int i = 0; i < 1; ++i)
    {
        double freq = 196.0 * pow (2, (7.0 * i) / 12.0);
        violinStrings.add (new ViolinString (freq, fs, 0, elastoPlastic));
    }
    
    setSize (appWidth, appHeight);
    
    Timer::startTimerHz (15);
    
    
    for (auto violinString : violinStrings)
    {
        addAndMakeVisible (violinString);
        std::cout << (violinString->getModel() == exponential ? "Exponential" : "ElastoPlastic") << std::endl;
    }
}

void MainComponent::getNextAudioBlock (const AudioSourceChannelInfo& bufferToFill)
{
    // Your audio-processing code goes here!

    // For more details, see the help for AudioProcessor::getNextAudioBlock()

    // Right now we are not producing any data, in which case we need to clear the buffer
    // (to prevent the output of random noise)
    float *const channelData1 = bufferToFill.buffer->getWritePointer(0, bufferToFill.startSample);
    float *const channelData2 = bufferToFill.buffer->getWritePointer(1, bufferToFill.startSample);
    
//    float output{0.0, 0.0};
    for (int i = 0; i < bufferSize; ++i)
    {
        float output = 0.0;
        for (auto violinString : violinStrings)
        {
            violinString->bow();
            output = output + violinString->getOutput(0.8) * (violinString->getModel() == exponential ? 800 : 800); // comment everything after "getOutput(0.25)" out when debugging on a sample-by-sample basis (vs. matlab)
            violinString->analyseEnergy();
            violinString->updateUVectors();
        }
        
        // Do this after uvectors have been updated
        if (updateGridFlag)
        {
            violinStrings[0]->updateGrid (waveSpeedVal);
            updateGridFlag = false;
        }
        channelData1[i] = clip(output);
        channelData2[i] = clip(output);
//        std::cout <<"Buffer sample " << i << ": " <<  output << std::endl;
    }
//    std::cout << violinStrings[0]->getzBA() << std::endl;
}

void MainComponent::releaseResources()
{
    // This will be called when the audio device stops, or when it is being
    // restarted due to a setting change.

    // For more details, see the help for AudioProcessor::releaseResources()
}

//==============================================================================
void MainComponent::paint (Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));
    
}

void MainComponent::resized()
{
    // This is called when the MainContentComponent is resized.
    // If you add any child components, this is where you should
    // update their positions.
    
    Rectangle<int> totalArea = getLocalBounds();

//    if (showControls)
//    {
        Rectangle<int> controlsRect = totalArea.removeFromRight(controlsWidth);
    
        waveSpeedLabel.setBounds (controlsRect.removeFromTop(30));
        waveSpeedSlider.setBounds (controlsRect.removeFromTop(controlsWidth * 0.8));
//
//        controlsRect.removeFromTop(20);
//        scaleLabel.setBounds (controlsRect.removeFromTop(30));
//        scaleGraphics.setBounds (controlsRect.removeFromTop (controlsWidth * 0.8));
//
//        controlsRect.removeFromTop(20);
//        noiseLabel.setBounds (controlsRect.removeFromTop(30));
//        noiseFactor.setBounds (controlsRect.removeFromTop (controlsWidth * 0.8));
//
//        overrideNoiseButton.setBounds (controlsRect.removeFromTop(50));
//        toggleGraphics.setBounds (controlsRect.removeFromTop(50));
//    }
    
    int i = 0;
    int div = appHeight / violinStrings.size();
    for (auto violinString : violinStrings)
    {
        violinString->setBounds (totalArea.removeFromTop(div));
        ++i;
    }
    
    int horSizeData = 300;

}

void MainComponent::timerCallback()
{
    repaint();
}

double MainComponent::clip (double output, double min, double max)
{
    if (output > max)
    {
        return output = max;
    }
    else if (output < min)
    {
        return output = min;
    }
    return output;
}

void MainComponent::sliderValueChanged (Slider* slider)
{
    if (slider == &waveSpeedSlider)
    {
        updateGridFlag = true;
        waveSpeedVal = slider->getValue();
    }
}

void MainComponent::buttonClicked (Button* button)
{
    
}
