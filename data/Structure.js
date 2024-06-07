class Command {
    constructor(title, isEnabled, serverCall) {
        this.title = title;
        this.isEnabled = isEnabled;
        this.serverCall = serverCall; // Empty serverCall implies a section header
    }

    setEnabled(isEnabled) {
        return new Command(this.title, isEnabled, this.serverCall);
    }
}

class Parameter {
    constructor(title, isSub, type, options, serverTag, value) {
        this.title = title;
        this.isSub = isSub;
        this.type = type;
        this.options = options;
        this.serverTag = serverTag;
        this.value = value;
    }

    setValue(newValue) {
        return new Parameter(this.title, this.isSub, this.type, this.options, this.serverTag, newValue);
    }
}

module.exports = {Command, Parameter};