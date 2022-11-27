import React, { Component } from "react";
import Form from "../forms/Form";
import AppfetchOb from "./API_act_fetch";
import AppApi from "./API_handler";

class AppLand extends Component {
  constructor(props) {
    super(props);
    this.state = {
      post: {
        name: ""
      },
      jobs: [],
      showForm: true
    };
  }


  //have 3 state variables, showForm will render the form, post(name) and jobs [] unknown

  handleChange = e => {
    const { name, value } = e.target;
    // console.log(name,value,'nameValue')
    this.setState(prevState => ({
      post: {[name]: value },
    }));
    
  };
  //When a value will be added to the user form, this change function will be called and then new value will be added to prevstate ('' in default firsttime), e.g. [name]: value is, store the value in it.
    /*
  Doc 16 Nov 22:
  
  state has three componenets
    post: {name: ''}
    jobs: [],
    showform: true

  if showform state 
    return 
    flex ob
        two flexes 50-50 width
            flex 1 
                <form>
                button example
            flex 2 
                <image carousal>
                inactive
  
  else
    it checks if 0th object of jobs has name vakue of string or integer type and coordingle hits the URL 
    


    how form works:

        form is passed handlechange, post and handlesubmit
            on every text written in field, it will add letter (handlechange)
            handlechange will add and remove text from the post.name state (concatenates to previous value it holds)

            on clicking handlesubmit, jobs list is created with single elem holding object with key as name and value as user string

  */

  handleSubmit = e => {
    e.preventDefault();
    this.setState(prevState => ({
      jobs: [prevState.post],
      showForm: false
    }));
  };

  example = e => {
    e.preventDefault();
    this.setState(prevState => ({
      jobs: [{name :"56417"}],
      showForm: false
    }));
  };

  render() {

    let exButton = (
      <div className="button">
        <button className="btn" onClick={this.example} style={{background:"#DCAE1d",color:"White"}}>
          Example
        </button>
      </div>
    )

    // console.log("JOBS LOG",this.state.jobs)
    // const urlfetch = "http://172.16.3.146/nextrap/";
    // const urlfetch = "http://localhost:8001";
    // const urlfetch = "http://14.139.227.206/nextrap";

    const urlfetch = "http://localhost:8000";
    // const urlfetch = "http://172.16.46.73/enactdb/"
    // console.log("landingpage", urlfetch, this.state);
    if (this.state.showForm) {
      return (
        <div className="d-flex h-100" style={{ alignItems:"center", verticalAlign: "middle"}}>
          {/* style={{alignSelf: "center", alignItems: "center", borderStyle: "double",  borderColor: "blue"}} */}
          <div className="d-flex w-50 p-0 m-0">
            {/* this is not working */}
            <div className="border">
              <Form
                handleChange={this.handleChange}
                post={this.state.post}
                handleSubmit={this.handleSubmit}
              />
              {exButton}
            </div>
          </div>
          <div className="w-50 h-100 p-0 m-0" style={{ borderRadius: "50% 0% 0% 50% / 50% 0% 0% 50%", 
        background: 'linear-gradient(to right, #E9E9E9, #E8E8E8'}}>
            <div className="position-relative p-5 m-5" align="center">
              <h3>will add image carousla here with hyperlink to about nomenclaure and turtrials</h3>
            </div>
          </div>
        </div>
      );
    } else {
      if (isNaN(this.state.jobs[0].name)) {
        // It will return true if not a number else false, means inside if its a string
        // console.log("insidestring funct", this.state.jobs[0], isNaN(this.state.jobs[0]));
        return (
          <div className="post-container">
            {this.state.jobs.map((job, index) => (
              <div key={index}>
                <AppApi
                  key={index}
                  urlprefix={urlfetch}
                  urlfetch={`${urlfetch}/name/${job.name}/`}
                />
              </div>
            ))}
          </div>
        );
      } else {
        // console.log("insidenumeric funct");
        return (
          <div className="post-container">
            {this.state.jobs.map((job, index) => (
              <div key={index}>
                <AppfetchOb
                  key={index}
                  urlfetch={`${urlfetch}/ncbid/${job.name}/`}
                />     
              </div>
            ))}
          </div>
        );
      }
    }
  }
}

// text input to the form is added to state "post.name" and on click 
// this information added to the name varible then will be pushed to the jobs stack and then the showForm state will be toggled false.

// if added identifier will be alphabetical, then 
    // isNaN(this.state.jobs[0].name) will be true and will go to first renedered block, and iwll try to fetch url with name job.name <component will be 'AppApi'
    // otherwse if numerical then the else block , delaing here can be tricky as if url isnt fetched and other details, this componenet may go awry. coponent will be 'AppfetchOb'

export default AppLand;
