#pragma once
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkObjectFactory.h>
#include <CShapeDisplay.h>

class KeyPressInteractorStyle : public vtkInteractorStyleTrackballCamera
{
  public:
    static KeyPressInteractorStyle* New();
    vtkTypeMacro(KeyPressInteractorStyle, vtkInteractorStyleTrackballCamera);

    CShapeDisplay * m_display;
    
    void SetShape(CShapeDisplay *temp)
    {
      m_display = temp;
    }

    
    virtual void OnKeyPress() 
    {
      // Get the keypress
      vtkRenderWindowInteractor *rwi = this->Interactor;
      std::string key = rwi->GetKeySym();
 
      // Output the key that was pressed
      std::cout << "Pressed " << key << std::endl;
 
      // Handle an arrow key
      if(key == "Up")
      {
        // std::cout << "The up arrow was pressed." << std::endl;
        m_display->key_up();
      }
 
      // Handle a "normal" key
      if(key == "Down")
      {
        // std::cout << "The a key was pressed." << std::endl;
        m_display->key_down();
      }

      if (key == "minus")
      {
        m_display->key_minus();
      }

      if (key == "equal")
      {
        m_display->key_equal();
      }

      if (key == "R" || key == "r")
      {
        m_display->key_r();
      }

      if (key == "f" || key == "F")
      {
        m_display->key_f();
      }

      if (key == "bracketleft")
      {
        m_display->key_bracketleft();
      }

      if (key == "bracketright")
      {
        m_display->key_bracketright();
      }

      if (key == "o")
      {
        m_display->key_o();
      }
      if (key=="p")
      {
        m_display->key_p();
      }
      
      if (key == "1" || key =="2" || key == "3" || key=="4"|| key== "5"||key=="0")
      {
        m_display->key_change_viewport(key);
      }
      // Forward events
      vtkInteractorStyleTrackballCamera::OnKeyPress();
    }
 
};











