class Parameter {
    constructor(title, type, options, placeholder, value) {
        this.title = title;
        this.type = type;
        this.options = options;
        this.placeholder = placeholder;
        this.value = value;
    }
}

class Executable {
    constructor(isEnabled, title, checkFile, environment, specialCase, parameters, command, downloadables) {
        this.isEnabled = isEnabled;
        this.title = title;
        this.checkFile = checkFile;
        this.environment = environment;
        this.specialCase = specialCase;
        this.parameters = parameters;
        this.command = command;
        this.downloadables = downloadables;
    }
}

class Downloadable {
    constructor(enabled, isVisual, label, path) {
        this.enabled = enabled;
        this.isVisual = isVisual;
        this.label = label;
        this.path = path;
    }
}

module.exports = {Parameter, Executable, Downloadable};