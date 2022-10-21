import React, { Component } from "react";
// import Headercon from "./components/header";
// import Footercon from "./components/footer";
import Navbar from "./components/navbar";
// import logo from "./nextrap.png";
//import AppLand from "./landing_api/Landing_page";

class PlaceHolder extends Component {
  
  render() {
    return (
      <div id="body-wrap">

        
        {/* This component is named navbar, however its not only renders navbar but also the regions thereafter */}
        <Navbar />

        {/* footer content below */}
      </div>
    );
  }
}

// const footer = {
//   background: "#ccff99",
//   textAlign: "center",
//   paddingTop: "2px",
//   paddingBottom: "2px",
//   color: "black",
//   font: "1em Arial, Tahoma, sans-serif",
//   overflow: "auto",
// }
export default PlaceHolder;
//style={{background:"linear-gradient(to bottom right, blue, green)"}}

//https://mdbootstrap.com/docs/react/extended/gradients/