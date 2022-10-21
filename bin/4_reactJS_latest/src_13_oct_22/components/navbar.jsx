import React, { Component } from "react";
import About from "./about";
import Tutorial from "./tutorial";
import AppLand from "../landing_api/Landing_page";
import Nomenclature from "./nomenclature";
//import Badge from 'react-bootstrap/Badge';
import HolderFront from '../gene_render/HolderFront';

class Navbar extends Component {
  constructor(props) {
    super(props);
    this.state = {
      home: true,
      about: false,
      tutorial: false,
      nomenclature: false
    };
  }
  toggle_home = e =>
    this.setState({
      about: false,
      tutorial: false,
      nomenclature: false,
      home: true
    });

  toggle_name = e =>
    this.setState({
      home: false,
      about: false,
      tutorial: false,
      nomenclature: true
    });

  toggle_about = e =>
    this.setState({
      home: false,
      tutorial: false,
      nomenclature: false,
      about: true
    });

  toggle_tutorial = e =>
    this.setState({
      home: false,
      about: false,
      nomenclature: false,
      tutorial: true
    });

  render() {

    return (
      <div className="container p-0">
        {/* style={{height: "840px"}} */}
        {/* style={{minHeight:"520px", maxHeight:"840px"}} */}
        <div className="container position-sticky fixed-top p-0">
          <nav className="navbar navbar-expand-md justify-content-center bg-dark"
          style={navGradient}
          > 
          {/*navbar navbar-expand-md justify-content-center m-0 navbar-dark bg-dark */}
            <div className="container-fluid">
              <div >
                  <a 
                  className="navbar-brand d-flex w-50 me-auto nav-link active" 
                  href="http://localhost:3000"
                  onClick={this.toggle_home}
                  style= {this.state.home ? activeHome : inactiveHome}
                  >ENACTdb</a>
              </div>
              <div className="navbar-collapse collapse w-100 mt-5 p-0">
                <ul className="nav navbar-nav ms-auto w-100 justify-content-end">
                    <li className="nav-item">
                      <a
                        className="nav-link"
                        href="#"
                        onClick={this.toggle_about}
                        style= {this.state.about ? active : inactive}
                      >About</a>
                    </li>
                    <li className="nav-item">
                    <a
                        className="nav-link"
                        href="#"
                        onClick={this.toggle_name}
                        style= {this.state.nomenclature ? active : inactive}
                      >Exon-nomenclature</a>
                    </li>
                    <li className="nav-item">
                    <a
                      className="nav-link"
                      href="#"
                      onClick={this.toggle_tutorial}
                      style= {this.state.tutorial ? active : inactive}
                    >Tutorial</a>
                    </li>
                </ul>
              </div>
            </div>
          </nav>
          </div>
        <div className="container">
          <div className="container py-3 border">
            {this.state.about ? <About/> : ' '}
            {this.state.nomenclature ? <Nomenclature/> : ' '}
            {this.state.home ? <AppLand/> : ' '}
            {/* {this.state.home ? <HolderFront/> : ' '} */}
            {this.state.tutorial ? <Tutorial/> : ' '}
            {console.log('InThere')}
          </div>
          </div>
        <div className="container p-0">
          <div className='bg-dark text-center text-lg-left'
              style={{ background: 'linear-gradient(to right, #2d545e, #12343b' }}
              >
                <div className='text-center p-3 text-light' style={{textIndent: "50px"}}>
                    Server developed and maintained by shashibp-lab: "
                    <a className='text-light' href='https://shashibp-lab.github.io/' >
                      https://shashibp-lab.github.io/"
                    </a>
                  {/* <a href="https://shashibp-lab.github.io" target="_blank"> shashibp-lab</a> */}
                  <br />
                  &copy; {new Date().getFullYear()} Copyright:{' '} Department of Biological Sciences, IISER Mohali, Sector 81, SAS
                  Nagar, 140306, Punjab, India
                </div>
          </div>
        </div>
      </div> 
    );
  }
}

const boldStyle = {
  color: "white",
  fontWeight: "thick"
};
const navGradient = {
  // background:"linear-gradient(to bottom right, #2d545e, #12343b)"
  background: "linear-gradient(to right, #52658f, #333a56",
  // style={{ background: 'linear-gradient(to right, #2d545e, #12343b' }} // blusih green 1st
  // style={{ background: 'linear-gradient(to right, rgba(102, 126, 234, 0.5), rgba(118, 75, 162, 0.5))' }}
  top: 0,
  // background: "transparent"
} 

const active = {
  color: "#A9A9A9",
  fontWeight: "thick",
  fontSize: 30,
  verticalAlign: "middle"
};
const inactive = {
  color: "white", //carbon
  fontWeight: "thick",
  fontSize: 30,
  verticalAlign: "middle"
};
const activeHome = {
  color: "white",
  fontFamily: "monoton",
  fontWeight: "thick",
  fontSize: 75,
  verticalAlign: "middle"
};
const inactiveHome = {
  color: "#A9A9A9", //carbon
  fontWeight: "thick",
  fontFamily: "monoton",
  fontSize: 30,
  verticalAlign: "middle"
};
const background = {
  //background: "#F7FSE6" // cant be distinguished from the normal white, no contrast actuially
  background: "#98878F",
  // background: "black"
  // minHeight: "600px",
  // maxHeight: "1260px",
  verticalAlign: "middle"

}
export default Navbar;
