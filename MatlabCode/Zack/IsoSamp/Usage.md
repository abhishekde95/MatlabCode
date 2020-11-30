# IsoSamp Background

The purpose of IsoSamp is to take a mathematical description (model) of a
subject's iso-detection performance, sample points on its surface, and present
those stimuli (the points) to an observer. The models are fit using the
subject's perception of a small set of stimuli. These models are then used to
"fill in the gaps," i.e., predict how well the subject can see the stimuli not
shown during the behavioral experiments. Once we have the model, we can present
the larger set of stimuli to neurons in rapid succession.

# IsoSamp Usage

## Before you start

You must run the Grating paradigm before using IsoSamp. This is because IsoSamp
inherits the preferred spatial frequency and orientation of the neuron from
Grating. This is crucial for DTNT-based stimuli, but not for LMTF stimuli. If
you are using IsoSamp to present LMTF stimuli, then you can skip running Grating
and hardcode the spatial frequency and orientation via (execute the following on
the Slave computer):

```matlab
>> global gl
>> gl.prefsf = 1;
>> gl.preforient = 0;
>> MasterSlaveProgram
```

## Online machine

**IsoSamp requires SortClient with data transfer started!**

Run `MasterOnlineProgram` like usual. When you run the paradigm for the first
time on REX or hit "Reset States," a file dialog will appear to select what kind
of stimuli to sample and send to REX (the DTNT- and LMTF-based stimuli are
currently supported).

The Online program will then prompt you to select which subject's parameters to
use. The subjects are identified by the first one or two letters of their name.

You will only need to select the module and subject once. There will be messages
that remind you whose model parameters are in use. If you need to change either
the module or the subject, `ctrl-c` the Online program, execute `clear all`,
restart `MasterOnlineProgram`, and click on "Reset States" within REX.

## REX

There are several menus in the IsoSamp REX paradigm containing variables that
change aspects of the experiment.

### Fixation and Gabor

These are self explanatory and not unique to this paradigm.

### Behavior

- There is a "behavior" flag that will present the stimulus either in the
receptive field or at the RF's position flipped about the y-axis. Then, the
subject reports which side the stimulus appeared with a saccade to one of two
targets. This is in contrast to the default mode of operation where the stimulus
always appears in the receptive field and the subject must maintain fixation.

- The "use previous stimuli" flag will signal the Online code to return the same
(L,M,S,TF) tuples as the previous run. This is useful if stimulus generation
involves some sort of randomization, but you wish to reproduce the exact
stimulus set. This flag is ignored if the receptive field changed, issuing a
warning on the Online side.

### Stimuli

- The target number of stims is a variable (**must be > 0 in order for the
experiment to begin**) that tells the Online code to attempt to return that
number of unique stimuli.
- Number of blocks: number of stimulus repeats taken in a block-wise fashion.
- Special case is a variable that signals the LMTF module to take specific
subsets of stimuli, ignoring the "Advanced" parameters below. They are:
  1. luminance;
  2. red-green;
  3. 1 and 2;
  4. 1 and 2 of both polarities.
- Blank trials flags whether a null stimulus should exist in the stimulus set.

### Advanced

The user can specify a subregion of the surface to pick stimuli from. This is
parameterized in cylindrical coordinates (r,theta,z):
- `r_lo <= 1 <= r_hi`, where 1 represents detection threshold in most cases,
defines an interval where an r is chosen uniformly.
- `0 <= theta_lo <= theta_hi < pi` defining a wedge in the LM-plane.
- `z_lo <= z_hi` restricts the z-coordinate (for example, z could be S cone
contrast or temporal frequency).

Note that the radius of each data point will be scaled by a scalar picked
uniformly in `[r_lo,r_hi]`.
