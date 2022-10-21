import React, { Component } from "react";
import PlaceHolder from "./placeholder";
import 'bootstrap/dist/css/bootstrap.min.css';
//import Track from "./exonVista/Track"
import ModuleExon from "./exonVista/ModuleExon";
  class App extends Component {
  
  render() {
    const reg=/H+/gd;
    const helixStr='-------------------HHHHHHHHHHHHHHHHHHHHHHHH----HH--H---HHH-'
    var result, indices = [];
    while ((result = reg.exec(helixStr)) ) {
      //indices.push(result.index);
      indices.push(result.indices);
  }
    //console.log(reg.exec(helixStr).indices[0]);
    //console.log('*',indices)
     {return <PlaceHolder />; }
    // { return (
    //   <> 
    //     { <ModuleExon/>  }
    //   </>)
    //  }
    //return <h2>"hello"</h2>
}
}
export default App;
